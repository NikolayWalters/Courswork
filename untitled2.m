clear all
tmax = 25;
step = 0.01;
tsol = 0:step:tmax;
ysol = zeros(1, length(tsol));
xsol = zeros(1, length(tsol));
errorx = zeros(1, length(tsol));
errory = zeros(1, length(tsol));
xsol(1) = 1;
ysol(1) = 1;
ysol1 = zeros(1, length(tsol));
xsol1 = zeros(1, length(tsol));
xsol1(1) = 1;
ysol1(1) = 1;
Aconst = 2;
Bconst = 6;
fxyt = @(tsol, xsol, ysol) Aconst - (Bconst * xsol) + xsol^2 * ysol - xsol;
gxyt = @(tsol, xsol, ysol) Bconst * xsol - xsol^2 * ysol;
for i = 1:(length(tsol)-1)
    k0 = fxyt(tsol(i), xsol(i), ysol(i));
    l0 = gxyt(tsol(i), xsol(i), ysol(i));
    k1 = fxyt(tsol(i)+0.5*step,xsol(i)+0.5*step*k0,ysol(i)+0.5*step*l0);
    l1 = gxyt(tsol(i)+0.5*step,xsol(i)+0.5*step*k0,ysol(i)+0.5*step*l0);
    k2 = fxyt((tsol(i)+0.5*step),(xsol(i)+step*(k0+k1)/4),(ysol(i)+step*(l0+l1)/4));
    l2 = gxyt((tsol(i)+0.5*step),(xsol(i)+step*(k0+k1)/4),(ysol(i)+step*(l0+l1)/4));
    k3 = fxyt((tsol(i)+step),(xsol(i)-step*(k1-2*k2)),(ysol(i)-step*(l1-2*l2)));
    l3 = gxyt((tsol(i)+step),(xsol(i)-step*(k1-2*k2)),(ysol(i)-step*(l1-2*l2)));
    k4 = fxyt((tsol(i)+(2*step/3)),(xsol(i)+step*(7*k0+10*k1+k3)/27),(ysol(i)+step*(7*l0+10*l1+l3)/27));
    l4 = gxyt((tsol(i)+(2*step/3)),(xsol(i)+step*(7*k0+10*k1+k3)/27),(ysol(i)+step*(7*l0+10*l1+l3)/27));
    k5 = fxyt((tsol(i)+0.2*step),(xsol(i) + step*(28*k0-125*k1+546*k2+54*k3-378*k4)/625),(ysol(i) + step*(28*l0-125*l1+546*l2+54*l3-378*k4)/625));
    l5 = gxyt((tsol(i)+0.2*step),(xsol(i) + step*(28*k0-125*k1+546*k2+54*k3-378*k4)/625),(ysol(i) + step*(28*l0-125*l1+546*l2+54*l3-378*k4)/625));
    xsol(i+1) = xsol(i) + (1/6)*step*(k0 + 4*k2 + k3);
    ysol(i+1) = ysol(i) + (1/6)*step*(l0 + 4*l2 + l3);
    xsol1(i+1) = xsol(i) + (1/336)*step*(14*k0+35*k3+162*k4+125*k5);
    ysol1(i+1) = ysol(i) + (1/336)*step*(14*l0+35*l3+162*l4+125*l5);
    
end
errorx = abs(xsol-xsol1);
errory = abs(ysol-ysol1);
%hold on
%plot(xsol, ysol,'DisplayName', 'Converged limit cycle');
%plot(xsol1, ysol1,'DisplayName', 'Accurate');
%title('RK4 approximation of the Brusselator in the phase plane (x,y)','FontSize',14);
%xlabel('x', 'FontSize',14);
%ylabel('y', 'FontSize',14);
%legend('show', 'FontSize',20);
%set(gca,'FontSize',20);
%hold off
hold on
plot(tsol, errorx,'DisplayName', 'X');
plot(tsol, errory,'DisplayName', 'y');
legend('show', 'FontSize',20);
