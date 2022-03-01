function lambda = calculateLambda(gamma,BC0,BC1,BC2)
M = BC0+gamma*BC1+gamma^2*BC2;
D = eig(M);
lambda = abs(max(real(D)));
end