clear all
tmax = 250;
step = 0.01;
tsol = 0:step:tmax;
ysol = zeros(1, length(tsol));
xsol = zeros(1, length(tsol));
xsol(1) = 0;
ysol(1) = 0;
Aconst = 1;
Bconst = 2.1;
fxyt = @(tsol, xsol, ysol) Aconst - (Bconst * xsol) + xsol^2 * ysol - xsol;
gxyt = @(tsol, xsol, ysol) Bconst * xsol - xsol^2 * ysol;
for i = 1:(length(tsol)-1)
    k0 = fxyt(tsol(i), xsol(i), ysol(i));
    l0 = gxyt(tsol(i), xsol(i), ysol(i));
    k1 = fxyt(tsol(i)+0.5*step,xsol(i)+0.5*step*k0,ysol(i)+0.5*step*l0);
    l1 = gxyt(tsol(i)+0.5*step,xsol(i)+0.5*step*k0,ysol(i)+0.5*step*l0);
    k2 = fxyt((tsol(i)+0.5*step),(xsol(i)+0.5*step*k1),(ysol(i)+0.5*step*l1));
    l2 = gxyt((tsol(i)+0.5*step),(xsol(i)+0.5*step*k1),(ysol(i)+0.5*step*l1));
    k3 = fxyt((tsol(i)+step),(xsol(i)+k2*step),(ysol(i)+l2*step));
    l3 = gxyt((tsol(i)+step),(xsol(i)+k2*step),(ysol(i)+l2*step));
    xsol(i+1) = xsol(i) + (1/6)*step*(k0 + 2*k1 + 2*k2 + k3);
    ysol(i+1) = ysol(i) + (1/6)*step*(l0 + 2*l1 + 2*l2 + l3);
end
hold on
plot(tsol, ysol,'DisplayName', 'A = 1, B = 2.5');
%plot(tsol, ysol,'r','DisplayName', 'Converged limit cycle');
%plot(xconst,yconst, 'DisplayName','Estimated boundary');
%plot(xcon,ycon,'DisplayName', 'Strict boundary');
title('RK4 y approximation of the Brusselator between t = 0 and 25','FontSize',14);
xlabel('t', 'FontSize',14);
ylabel('Approximation', 'FontSize',14);
legend('show', 'FontSize',20);
set(gca,'FontSize',20);

