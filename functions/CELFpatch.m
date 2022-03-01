function [res_new] = CELFpatch(res,params,b1plus,mask,varargin)
%CELFpatch estimates T1, T2, off-resonance and banding free image with
%ellipse fitting constrained center on a line for a 3x3 window
%   ----inputs----
%       res:     results obtained from voxel-wise CELF
%       params:  struct including sequence parameters
%       b1plus:  B1+ map
%       mask:    mask 
% 
%   ----output----
%       res_new: struct including refined estimation results
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
maskarr = reshape(mask,[N,1]); % vectorize the mask
b1arr = reshape(b1plus,[N,1]); % vectorize the B1+ map

if(size(b1arr,2) ~= 1)
    b1arr = reshape(b1plus,[N,1]);
end

res_new = res;

%% Iterate over all voxels that raised a symmetric check flag and inside the mask
for i = 2:Nx-1
    for j = 2:Ny-1
        if(res.symStatus(i,j))
            if(~any(mask(i-1:i+1,j-1:j+1) == 0))
                
                ind = sub2ind([Nx,Ny],i,j);
                
                %% Aggregate all data points in a 3x3 window
                
                [ii,jj] = meshgrid(i-1:i+1,j-1:j+1);
                meshind = sub2ind([Nx,Ny],ii,jj);
                indlist = [];
                indlist2 = [];
                for k = 1:9
                    if(maskarr(meshind(k))==0)
                        indlist = [indlist, (k-1)*Np+1:k*Np];
                        indlist2 = [indlist2, k];
                    end
                end
                
                a = res.rotS(i-1:i+1,j-1:j+1,:);
                md = [real(a(:)),imag(a(:))];
                md(indlist,:) = [];
                
                allcpl = res.Meff(i-1:i+1,j-1:j+1);
                tempcpl = allcpl(:);
                tempcpl(indlist2,:) = [];
                cpl = mean(tempcpl);
                
                %% Constrained ellipse fitting with center on a line
                
                % coefficients of the rotated ellipse
                rotC = ellipseFitConstLs(md,norm(cpl));
                res_new.rotCarr(ind,:) = rotC;
                
                % center and semi-axes of the rotated ellipse
                [cent,A,B] = findEllipseCenterAndAxes(rotC);
                
                %% Analytical solutions for parameters T1 and T2
                
                % find analytic solution for parameters a and b
                [res_new.a(ind),res_new.b(ind),~] = findAnalyticSol(norm(cent),norm(A),norm(B));
                
                % estimate T1 and T2
                [res_new.T1arr(ind),res_new.T2arr(ind)] = estimateT1T2(res_new.a(ind),res_new.b(ind),TR,FArad*b1arr(ind));
                
                %% Estimation of the local off-resonance delta_f0
                
                costn = (md(:,1)-norm(cent))/norm(A);
                costhn = (costn-res_new.b(ind))./(res_new.b(ind)*costn-1);
                
                cosdthn = cos(repmat(deltaThetas,[1,numel(tempcpl)]));
                sindthn = sin(repmat(deltaThetas,[1,numel(tempcpl)]));
                
                K = [cosdthn',sindthn']\costhn;
                
                theta0 = atan2(K(2),K(1));
                
                res_new.theta0arr(ind) = theta0;
                
                % off-resonance estimation
                res_new.f0arr(ind) = theta0/(2*pi*TR)*10^3;
                
                %% Estimate banding-free bSSFP signal
                
                res_new.IBFarr(ind) = res_new.Meffarr(ind).*(1+res_new.a(ind))./(1+res_new.b(ind));

            end
        end
    end
end

IBF = reshape(res_new.IBFarr,[Nx, Ny]);
T1 = reshape(res_new.T1arr,[Nx, Ny]);
T2 = reshape(res_new.T2arr,[Nx, Ny]);
f0 = reshape(res_new.f0arr,[Nx, Ny]);
Meff = reshape(res_new.Meffarr,[Nx, Ny]);

res_new.T1 = T1;
res_new.T2 = T2;
res_new.f0 = f0;
res_new.IBF = IBF;
res_new.MEff = Meff;

end