function MinObj = P4(x,pf)
global ans
sz = size(x);
num = sz(1,2);
term1 = 0;
term2 = 1;
term3 = 0;
for i = 1:num
    term1 = term1 + (cos(x(1,i)))^4;
    term2 = (term2*((cos(x(1,i)))^2));
    term3 = term3 + i*((x(1,i))^2);
end
term4 = (abs(term1-2*term2));
term5 = sqrt(term3);
ans = -term4/term5;
sc1 = 1;
sc2 = 0;
for i = 1:num
    sc1 = sc1*x(1,i);
    sc2 = sc2 + x(1,i);
end
s1 = 0.75-sc1;
s2 = sc2-((15*num)/2);
if ans == -Inf
    pf = Inf;
end
pf1 = pf*10^12;
pf2 = pf*10^6;
MinObj = (ans) + pf1*((max(0,s1))^2) + pf2*((max(0,s2))^2);
