function MinObj = P3(x,pf)
global ans
sz = size(x);
ans = (x(1,1))^2 + (0.5*x(1,1)) + (3*x(1,1)*x(1,2)) + (5*(x(1,2))^2);   
s1 = 3*x(1,1) + 2*x(1,2) + 2;
s2 = 15*x(1,1) - 3*x(1,2) - 1;
pf1 = pf*10^8;
pf2 = pf*10;
MinObj = (ans) + pf1*((max(0,s1))^2) + pf2*((max(0,s2))^2);
