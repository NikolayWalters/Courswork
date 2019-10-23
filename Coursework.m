clear all
tmax = 25;  % t range
step = 0.01; % h value
tsol = 0:step:tmax; % t array
ysol = zeros(1, length(tsol)); % y array
xsol = zeros(1, length(tsol)); % x array
ysol(1) = 11; % initial y
xsol(1) = 11; % initial x
Aconst = 2; % A constant
Bconst = 6; % B constant
fxyt = @(tsol, xsol, ysol) Aconst - (Bconst * xsol) + xsol^2 * ysol - xsol; % f(t,x,y)
gxyt = @(tsol, xsol, ysol) Bconst * xsol - xsol^2 * ysol; % g(t,x,y)
for i = 1:(length(tsol)-1) % RK4
    k0 = fxyt(tsol(i), xsol(i), ysol(i));
    l0 = gxyt(tsol(i), xsol(i), ysol(i));
    k1 = fxyt(tsol(i)+0.5*step,xsol(i)+0.5*step*k0,ysol(i)+0.5*step*l0);
    l1 = gxyt(tsol(i)+0.5*step,xsol(i)+0.5*step*k0,ysol(i)+0.5*step*l0);
    k2 = fxyt((tsol(i)+0.5*step),(xsol(i)+0.5*step*k1),(ysol(i)+0.5*step*l1));
    l2 = gxyt((tsol(i)+0.5*step),(xsol(i)+0.5*step*k1),(ysol(i)+0.5*step*l1));
    k3 = fxyt((tsol(i)+step),(xsol(i)+k2*step),(ysol(i)+l2*step));
    l3 = gxyt((tsol(i)+step),(xsol(i)+k2*step),(ysol(i)+l2*step));
    xsol(i+1) = xsol(i) + (1/6)*step*(k0 + 2*k1 + 2*k2 + k3); % x_{n+1}
    ysol(i+1) = ysol(i) + (1/6)*step*(l0 + 2*l1 + 2*l2 + l3); % y_{n+1}
end
%hold on
%title('RK4 approximation of the Brusselator between t = 0 and 5','FontSize',14);
%xlabel('t', 'FontSize',14);
%ylabel('Approximation', 'FontSize',14);
%plot(tsol, xsol,'DisplayName', 'x');
%plot(tsol, ysol, 'DisplayName', 'y');
%set(gca,'FontSize',20);
%legend('show', 'FontSize',20);
%hold off
hold on
plot(xsol, ysol, 'DisplayName', 'x_0 = 11, y_0 = 11');
title('RK4 approximation of the Brusselator in the phase plane (x,y)','FontSize',14);
xlabel('x', 'FontSize',14);
ylabel('y', 'FontSize',14);
legend('show', 'FontSize',20);
set(gca,'FontSize',20);