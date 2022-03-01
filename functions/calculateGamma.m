function [gamma,deltaFlag] = calculateGamma(C0,C1,C2,B)
t1 = C0(2,1)*B(1,2)+C0(1,2)*B(2,1);
t2 = C1(2,1)*B(1,2)+C1(1,2)*B(2,1);
t3 = C2(1,1)*C0(2,2)-C1(1,2)*C1(2,1);
t4 = C1(1,1)*C0(2,2)-C1(2,1)*C0(1,2)-C1(1,2)*C0(2,1);
t5 = C0(1,1)*C0(2,2)-C0(1,2)*C0(2,1);

t1 = t1*1e6;
t2 = t2*1e6;
t3 = t3*1e6;
t4 = t4*1e6;
t5 = t5*1e6;

delta = (t1*t2-8*t4)^2-4*(t2^2-16*t3)*(t1*t2*t4-t2^2*t5-4*t4^2)/t3;

if(delta<0)
    gam(1) = (-(t1*t2-8*t4) + sqrt(delta))/(t2^2-16*t3);
    gam(2) = (-(t1*t2-8*t4) - sqrt(delta))/(t2^2-16*t3);
    gam = abs(gam);
    deltaFlag = 1;
else
    gam(1) = (-(t1*t2-8*t4) + sqrt(delta))/(t2^2-16*t3);
    gam(2) = (-(t1*t2-8*t4) - sqrt(delta))/(t2^2-16*t3);
    deltaFlag = 0;
end
gamma = gam;
end