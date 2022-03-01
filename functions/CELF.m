function [results] = CELF(S,params,b1plus,mask,varargin)
%CELF estimates T1, T2, off-resonance and banding free image with
%ellipse fitting constrained center on a line
%   ----inputs----
%       S:       bSSFP image -complex valued, Np phase cycled, size of (Nx,Ny,Np)
%       params:  struct including sequence parameters
%       b1plus:  B1+ map
%       mask:    mask 
% 
%   ----output----
%       results: struct including estimation results
%--------------------------------------------------------------------------

%% Variables
if nargin ~= 4
    disp('Wrong number of inputs!');
    return;
end

TR = params.TR;
FArad = params.FArad;
deltaThetas = params.deltaThetas;
Nx = params.Nx;
Ny = params.Ny;
Np = params.Np;

N = Nx*Ny; % total number of voxels 
Sarr = reshape(S,[N,Np]); % vectorize the data
maskarr = reshape(mask,[N,1]); % vectorize the mask
b1arr = reshape(b1plus,[N,1]); % vectorize the B1+ map

a = zeros(N,1);
b = zeros(N,1);
T1arr = zeros(N,1);
T2arr = zeros(N,1);
f0arr = zeros(N,1);
theta0arr = zeros(N,1);
IBFarr = zeros(N,1);
rotCarr = zeros(N,6);
Meffarr = zeros(N,1);
symStatusarr = zeros(N,1);
rotSarr = zeros(N,Np);

%% Iterate over all voxels inside the mask
for i = 1:N
    
    if( maskarr(i) )
        
        % Split data into real and imaginary parts
        data = [(real(Sarr(i,:)) ); imag(Sarr(i,:))]';
        
        %% Step1: Find cross point passing through central line of the ellipse via ESM
        
        cpl = findCenterLineESM(data);
        Meffarr(i) = norm(cpl);
        
        %% Step2: Rotation of the ellipse to initial vertical conic form
        
        % rotation angle of the ellipse
        rotAngle = atan2(cpl(2),cpl(1));
        
        % rotated data points
        rotData = rotateDataPoints(data,-rotAngle);
        rotSarr(i,:) = rotData(:,1)' + 1i*rotData(:,2)';
        
        % Check if there are symmetric data points
        if(size(rotData,1)<6)
            deg = 12;
            symStatusarr(i) = symmetryCheck(rotData,norm(cent),deg);
        end
        
        %% Step3: Constrained ellipse fitting with center on a line
        
        % coefficients of the rotated ellipse
        rotC = ellipseFitConstLs(rotData,norm(cpl));
        rotCarr(i,:) = rotC;
        
        % center and semi-axes of the rotated ellipse
        [cent,A,B] = findEllipseCenterAndAxes(rotC);
        
        %% Step4: Analytical solutions for parameters T1 and T2
        
        % find analytic solution for parameters a and b
        [a(i),b(i),~] = findAnalyticSol(norm(cent),norm(A),norm(B));
        
        % estimate T1 and T2
        [T1arr(i),T2arr(i)] = estimateT1T2(a(i),b(i),TR,FArad*b1arr(i));
        
        %% Step5: Estimation of the local off-resonance delta_f0
        
        costn = (rotData(:,1)-norm(cent))/norm(A);
        costhn = (costn-b(i))./(b(i)*costn-1);
        
        cosdthn = cos(deltaThetas);
        sindthn = sin(deltaThetas);
        
        K = [cosdthn',sindthn']\costhn;
        
        theta0arr(i) = atan2(K(2),K(1));
        
        % off-resonance estimation
        f0arr(i) = (theta0arr(i))/(2*pi*TR)*10^3;
        
        %% Step6: Estimate banding-free bSSFP signal
        
        IBFarr(i) = Meffarr(i).*(1+a(i))./(1+b(i));
        
    else
        IBFarr(i) = 0;
        T1arr(i) = 0;
        T2arr(i) = 0;
        f0arr(i) = 0;
    end
end

IBF = reshape(IBFarr,[Nx, Ny]);
T1 = reshape(T1arr,[Nx, Ny]);
T2 = reshape(T2arr,[Nx, Ny]);
f0 = reshape(f0arr,[Nx, Ny]);
Meff = reshape(Meffarr,[Nx, Ny]);
rotS = reshape(rotSarr,[Nx,Ny,Np]);
symStatus = reshape(symStatusarr,[Nx, Ny]);

IBF(isinf(IBF)) = 0;
IBF(isnan(IBF)) = 0;
T1(isinf(T1)) = 0;
T1(isnan(T1)) = 0;
T2(isinf(T2)) = 0;
T2(isnan(T2)) = 0;
f0(isinf(f0)) = 0;
f0(isnan(f0)) = 0;
Meff(isinf(Meff)) = 0;
Meff(isnan(Meff)) = 0;

results.T1 = T1;
results.T2 = T2;
results.f0 = f0;
results.IBF = IBF;
results.MEff = Meff;

results.a = a;
results.b = b;
results.T1arr = T1arr;
results.T2arr = T2arr;
results.f0arr = f0arr;
results.theta0arr = theta0arr;
results.IBFarr = IBFarr;
results.rotCarr = rotCarr;
results.Meffarr = Meffarr;
results.symStatus = symStatus;
results.rotS = rotS;
end
