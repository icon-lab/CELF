function [a,b,Meff] = findAnalyticSol(xc,A,B)
%findAnalyticSol find analytic solution for parameters a,b,Meff
%   ---input--- 
%           xc: geometric center of the ellipse
%           A: semi-axes of the ellipse
%           B: semi-axes of the ellipse
%   ---output---
%           a: a coefficient of the complex signal
%           b: a coefficient of the complex signal
%           Meff: effective magnetization
%
%   where complex signal is modelled as:
%   I = Meff*(1-a*exp(i*theta))/(1-b*cos(theta))*exp(i*phase);
%   For detailed explanation please see appendix of following paper:
%   http://onlinelibrary.wiley.com/doi/10.1002/mrm.26717/pdf

b1 = (-xc*A + sqrt((xc*A).^2-(xc.^2+B^2)*(A^2-B^2)))./(xc.^2+B^2);
a1 = B/(xc.*sqrt(1-b1.^2) + b1*B);
Meff1 = xc.*(1-b1^2)/(1-a1*b1);

b2 = (xc*A + sqrt((xc*A).^2-(xc.^2+B^2)*(A^2-B^2)))./(xc.^2+B^2);
a2 = B/(xc.*sqrt(1-b2.^2) + b2*B);
Meff2 = xc.*(1-b2^2)/(1-a2*b2);

a = a1;
b = b1;
Meff = Meff1;

end

