function [C] = ellipseFitConstLs(data,cp)
%ELLIPSEFITCONSTLS find the coefficients of the ellipse via Constrained 
%Ellipse Fitting Formulation
%   ----inputs----
%       data: contains data point
%             Sx = data(:,1);
%             Sy = data(:,2);
%             data pairs should have -pi- degree between them and
%             should be ordered (eg.: 0-pi/4-pi/2-pi-5pi/4-3pi/2 for N=8)
%       cp:   Cross-point of the ellipse
% 
%   ----output----
%       C:    Coefficients of the ellipse          
%--------------------------------------------------------------------------

if all(data(:) == 0)
    C = zeros(6,1);
    return;
end

Sx = data(:,1);
Sy = data(:,2);

gammaBound = [0.5 1];

[x_c,u_opt,g] = constElFit(Sx,Sy,cp,gammaBound);

if(any(isnan(u_opt(:))))
    C = zeros(6,1);
    return;
end

C(1) = u_opt(1);
C(2) = 0;
C(3) = u_opt(2);
C(4) = -2*C(1)*x_c(1)-C(2)*x_c(2);
C(5) = -2*C(3)*x_c(2)-C(2)*x_c(1);
C(6) = C(1)*x_c(1)^2+C(3)*x_c(2)^2+C(2)*x_c(1)*x_c(2)+g(1);

end