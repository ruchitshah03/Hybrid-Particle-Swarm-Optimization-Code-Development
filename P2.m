function MinObj = P2(x)
global ans
sz = size(x);
num = sz(1,2);
ans = 0;
for i=1:num
    ans = ans + ((x(1,i))^2 - (10*cos(2*pi*x(1,i))));
end
MinObj = 10*num + ans;