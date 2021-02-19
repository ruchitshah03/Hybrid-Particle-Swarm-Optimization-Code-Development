function [c,ceq] = constraint(x)
sz = size(x);
num = sz(1,2);
sc1 = 1;
sc2 = 0;
for i = 1:num
    sc1 = sc1*x(1,i);
    sc2 = sc2 + x(1,i);
end
s1 = 0.75-sc1;
ceq = [];
s2 = sc2-((15*num)/2);
c = [s1, s2];
