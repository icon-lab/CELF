function [ cent,A,B ] = findEllipseCenterAndAxes(C)
%findEllipseCenterAndAxes find geometric center and semi-axes 
%of an ellipse
%   ---input--- 
%           C: coefficients of quadratic form of the ellipse
%   ---output---
%           xc: geometric center of the ellipse
%           A: semi-minor axes of the ellipse (vector)
%           B: semi-major axes of the ellipse (vector)
%   For detailed explanation please see link below:
%   https://www.geometrictools.com/Documentation/InformationAboutEllipses.pdf

a11 = C(1);
a12 = C(2)/2;
a22 = C(3);
b1 = C(4);
b2 = C(5);
c = C(6);

% center K=(k1,k2)
K = [a22*b1-a12*b2,a11*b2-a12*b1]./(2*(a12^2-a11*a22));
k1 = K(1);
k2 = K(2);

mu = 1/(a11*k1^2+2*a12*k1*k2+a22*k2^2-c);
m11 = mu*a11;
m12 = mu*a12;
m22 = mu*a22;

l1 = ((m11+m22)+sqrt((m11-m22)^2+4*m12^2))/2;
l2 = ((m11+m22)-sqrt((m11-m22)^2+4*m12^2))/2;

smin = 1/sqrt(l1); %semi-minor axis
smaj = 1/sqrt(l2); %semi-major axis

%minor axis direction
if(m11 >= m22)
    U1 = [l1-m22,m12];
    U1 = U1./norm(U1);
else
    U1 = [m12,l1-m11];
    U1 = U1./norm(U1);
end

%major axis direction
U2 = [-U1(2),U1(1)];

U1 = U1./norm(U1);
U2 = U2./norm(U2);

cent = K;
A = smin*U1;
B = smaj*U2;
end

