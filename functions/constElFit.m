function [x_c,u_opt,g] = constElFit(x,y,cp,gammaBound)
%CONSTELFIT performs constrained ellipse fitting to find the gamma given
%the data points and cross-point of the ellipse
%   ----inputs----
%       x, y:       Data points coordinates vectors
%       cp:         Cross-point of the ellipse
%       gammaBound: Vector for lower and upper bounds of gamma
% 
%   ----output----
%       x_c:        Center of the ellipse
%       u_opt:      Optimum parameter vector for ellipse parameters
%       g:          g
%--------------------------------------------------------------------------

% Some parts of this function are modified from the following reference:
% Waibel et al., Constrained Ellipse Fitting with Center on a Line,
% Journal of Mathematical Imaging and Vision, 2015.

%% Data Check
N = length(x);
T = [x.^2 x.*y y.^2 x y ones(N,1)];
if (rank(T,1e-8)==3 || cp==0)
    % Data points on a line
    x_c = [0 0];
    u_opt = [0 0];
    g = 0;
    return
end

%% Preparation 

B = [0 2; 2 0];
Z = eye(N)-ones(N)/N;

% centralized D0 and D1
D0 = [(x-cp).^2, y.^2];
D1 = [2*cp*(x-cp), y*0];

C0 = D0'*Z*D0;
C1 = -D0'*Z*D1-D1'*Z*D0;
C2 = D1'*Z*D1;

BC0 = B\C0;
BC1 = B\C1;
BC2 = B\C2;

gammaBound = gammaBound-1; % because of a centralized data

% Calculate gamma analytically
[gam,deltaFlag] = calculateGamma(C0,C1,C2,B);

% If analytical gamma is valid
if(~deltaFlag)
    lam = zeros(size(gam));
    for i=1:length(gam)
        lam(i) = calculateLambda(gam(i),BC0,BC1,BC2);
    end
    [~,idx]= min(lam);

    gamma = gam(idx);

    if(gamma<(gammaBound(1)))
       gamma = gammaBound(1);
    elseif(gamma>(gammaBound(2)))
       gamma = gammaBound(2);
    end
% Otherwise perform one dimensional search (modified from Waibel et al.)
else
    cnt = 0;
    maxIter = 100;
    isGlobalMinimumFound = false;

    while (~isGlobalMinimumFound && cnt<maxIter)

        cnt = cnt+1;

        [gamma,lambda,~,~] = fminbnd(@(gm) calculateLambda(gm,BC0,BC1,BC2),gammaBound(1),gammaBound(2),optimset('TolX',1e-8));
        
        %% Tunneling

        % Elimination of infinity-Eigenvalues
        P = [1 0; 0 cp^2];
        C2t = P'*C2*P;
        C1t = P'*C1*P;
        C0t = P'*C0*P;
        Bt = P'*B*P;
        C0t = C0t - lambda*Bt;

        C2b = C2t(1,1)-C1t(1,2)*C1t(2,1)/C0t(2,2);
        C1b = C1t(1,1)-C1t(1,2)*C0t(2,1)/C0t(2,2)-C0t(1,2)*C1t(2,1)/C0t(2,2);
        C0b = C0t(1,1)-C0t(1,2)*C0t(2,1)/C0t(2,2);
        rootsVec = eig([0 1; C0b C1b], [1 0;0 -C2b]);

        realRootsIdx = find(abs(imag(rootsVec)) < epsRootsImag & abs(rootsVec-gamma) > epsRootsMin);

        %% Treatment of Boundaries

        if numel(realRootsIdx) == 0
            isGlobalMinimumFound = true;

        elseif numel(realRootsIdx) == 1
            % Solution on given bounds
            break;

        elseif numel(realRootsIdx) == 3
            realRoots = rootsVec(realRootsIdx);
            rootsInsideIdx = find((realRoots > gammaBound(1)) & (realRoots < gammaBound(2)));
            if numel(rootsInsideIdx) == 1
                disp('Not possible')
            elseif numel(rootsInsideIdx) == 2
                gammaBoundNew(1) = min(real(realRoots(rootsInsideIdx)));
                gammaBoundNew(2) = max(real(realRoots(rootsInsideIdx)));
                disp('First solution on bound. New search bounds found inside.')
            end

        elseif numel(realRootsIdx) == 2
            gammaBoundNew(1) = min(real(rootsVec(realRootsIdx)));
            gammaBoundNew(2) = max(real(rootsVec(realRootsIdx)));

            if (abs(gammaBoundNew(2)-gammaBoundNew(1)) < epsRootsMin)

                lambdaNew = calculateLambda(gamma,BC0,BC1,BC2);

                if abs(lambdaNew - lambda) < epsFmin
                    disp(['Second global minimum found at gamma = ',num2str(gammaBoundNew(1))])
                else %Tunneling would lead to "perfectly fitting" non-ellipse 
                end
                break;
            else

                %Check if new bounds after tunneling are within given bounds
                gammaBoundCheck(1) = max(gammaBound(1),gammaBoundNew(1));
                gammaBoundCheck(2) = min(gammaBound(2),gammaBoundNew(2));

                if (gammaBoundCheck(1) > gammaBoundCheck(2))
                    disp('New Bounds outside given Bounds')
                    break;
                else
                    gammaBound = gammaBoundCheck;
                end
            end
        end
    end
end

%% Ellipse Parameters

% Find eigenvector u_opt
M = BC0+gamma*BC1+gamma^2*BC2;
[V,D] = eig(M);
idx = real(diag(D))==max(real(diag(D)));
u = real(V(:,idx));
u_opt = sign(u(1))*u/sqrt(u'*B*u);

% h_opt
h_opt = -ones(1,N)*(D0-gamma*D1)*u_opt/N;

% Find the center
x_c = cp*(gamma+1);

% g
g = h_opt - gamma^2*cp^2*u_opt(1);

end


