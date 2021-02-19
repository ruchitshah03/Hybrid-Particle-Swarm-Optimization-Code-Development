function [ans] = P5(x,pf)

m = 1;
k = x(1);
c = x(2);
omega = 0.1;
F0 = 1;

C = sqrt( (k-m*omega^2)^2 + (c*omega)^2 );
alpha = atan( c*omega/(k -m*omega^2));
omegad = sqrt( k/m - (c/(2*m))^2 );
A = - (F0/C) * cos(alpha);
B = -(F0/(C*omegad)) * ( omega*sin(alpha) + c/(2*m) * cos(alpha));

%%
% Calculate the analytical solution in the time interval [0,10] seconds at 500 equally space time intervals

f = fopen('MeasuredResponse.dat');
data = textscan(f, '%f %f');
data = cell2mat(data); % convert to matrix from cell array
fclose(f);

t = 0;
for i = 1:500
    term1 = A*cos(omegad*t) + B*sin(omegad*t);
    term2 = exp( - c*t/(2*m) );
    term3 = (F0/C)*cos(omega*t - alpha);
    u(i) = term1*term2 + term3;
    time(i) = t;
    t = t + 10/499;
    udot(i) = term2*(B*omegad*cos(omegad*t) - A*omegad*sin(omegad*t)) + (A*cos(omegad*t) + B*sin(omegad*t))*term2*(-c/(2*m)) - (F0/C)*omega*sin(omega*t-alpha);
    udotdot(i) = -term2*omegad^2*term1 - (c/(2*m))*term2*B*omegad*cos(omegad*t) + term2*(c/(2*m))*A*omegad*sin(omegad*t) + ((c/(2*m))^2)*term1*term2 - (c/(2*m))*term2*omegad*(B*cos(omegad*t)-A*sin(omegad*t)) - (F0/C)*(omega^2)*cos(omega*t-alpha);
    s1(i) = m*udotdot(i) + c*udot(i) + k*u(i) - F0*cos(omega*t);
    ans1(i) = (abs((u(i)'-data(i,2))));
    ans1(i) = ans1(i).^2;
end
pf1 = pf*1.1;
ans = sum(ans1) + pf1*(sum(s1.^2));
