function status = symmetryCheck(rotData,cent,deg)
N = size(rotData,1);
ang = zeros(N,1);
for i=1:N
    p1 = rotData(i,:);
    p2 = rotData(mod(i,N)+1,:);
    p3 = (p1+p2)/2;
    ang(i) = atan(p3(2)/(p3(1)-cent))*180/pi;
end
if(sum(abs(ang)<=deg)>=N/2)
    status = 1;
else
    status = 0;
end
end