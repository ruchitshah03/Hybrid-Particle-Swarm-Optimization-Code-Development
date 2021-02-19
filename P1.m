function MinObj = P1(x)
global ans
sz = size(x);
num = sz(1,2)-1;
ans = 0;
for i=1:num
    ans = ans + ((100*((x(1,i+1)-((x(1,i))^2))^2)) + ((1-(x(1,i)))^2));
end
MinObj = ans;