function [cpl] = findCenterLineESM(data)
%FINDCENTERLINEESM finds center line of the ellipse based on ESM
%   ---input--- 
%       data: contains data point
%             Sx = data(:,1);
%             Sy = data(:,2);
%             data pairs should have -pi- degree between them and
%             should be ordered (eg.: 0-pi/4-pi/2-pi-5pi/4-3pi/2 for N=8)
%   ---output---
%       cpl:  cross point of the ellipse - a point on the central line

Sx = data(:,1);
Sy = data(:,2);

Np = length(Sx); % number of phase cycles
NP = Np/2; % number of data pairs

Al = zeros(NP,2);
bl = zeros(NP,1);

for p = 1:NP
    x1 = Sx(p);
    x2 = Sx(p+NP);
    y1 = Sy(p);
    y2 = Sy(p+NP);
    Al(p,:) = [(y2-y1) , (x1-x2)];
    bl(p,1) = x1*(y2-y1)-y1*(x2-x1);
end

% cpl = Al\bl;
cpl = lsqminnorm(Al,bl);

end