function [ T1,T2 ] = estimateT1T2( a,b,TR,FA )
%estimateT1T2 estimates T1 and T2 from parameters a and b
%   ---input--- 
%           a: a coefficient of the complex signal
%           b: a coefficient of the complex signal
%           TR: repetition time (ms)
%           FA: flip angle (rad)
%   ---output---
%           T1: T1 time (ms)
%           T2: T2 time (ms)

cfa = cos(FA);

T1 = -TR/(log((a*(1+cfa-a*b*cfa)-b)/(a*(1+cfa-a*b)-b*cfa)));
T2 = -TR/log(a);

end

