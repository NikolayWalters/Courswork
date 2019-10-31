tmax = 25;
% a loop of h values would do much better here
% method below is equivalent to the code given in the report
step = 0.0001;
tsol = 0:step:tmax;
ysol = zeros(1, length(tsol));
xsol = zeros(1, length(tsol));
xsol(1) = 1;
ysol(1) = 1;
Aconst = 1.5;
Bconst = 2;
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
step = 0.001;
tsol0 = 0:step:tmax;
ysol0 = zeros(1, length(tsol0));
xsol0 = zeros(1, length(tsol0));
xsol0(1) = 1;
ysol0(1) = 1;
Aconst = 1.5;
Bconst = 2;
fxyt = @(tsol0, xsol0, ysol0) Aconst - (Bconst * xsol0) + xsol0^2 * ysol0 - xsol0;
gxyt = @(tsol0, xsol0, ysol0) Bconst * xsol0 - xsol0^2 * ysol0;
for i = 1:(length(tsol0)-1)
    k0 = fxyt(tsol0(i), xsol0(i), ysol0(i));
    l0 = gxyt(tsol0(i), xsol0(i), ysol0(i));
    k1 = fxyt(tsol0(i)+0.5*step,xsol0(i)+0.5*step*k0,ysol0(i)+0.5*step*l0);
    l1 = gxyt(tsol0(i)+0.5*step,xsol0(i)+0.5*step*k0,ysol0(i)+0.5*step*l0);
    k2 = fxyt((tsol0(i)+0.5*step),(xsol0(i)+0.5*step*k1),(ysol0(i)+0.5*step*l1));
    l2 = gxyt((tsol0(i)+0.5*step),(xsol0(i)+0.5*step*k1),(ysol0(i)+0.5*step*l1));
    k3 = fxyt((tsol0(i)+step),(xsol0(i)+k2*step),(ysol0(i)+l2*step));
    l3 = gxyt((tsol0(i)+step),(xsol0(i)+k2*step),(ysol0(i)+l2*step));
    xsol0(i+1) = xsol0(i) + (1/6)*step*(k0 + 2*k1 + 2*k2 + k3);
    ysol0(i+1) = ysol0(i) + (1/6)*step*(l0 + 2*l1 + 2*l2 + l3);
end
step = 0.01;
tsol1 = 0:step:tmax;
ysol1 = zeros(1, length(tsol1));
xsol1 = zeros(1, length(tsol1));
xsol1(1) = 1;
ysol1(1) = 1;
Aconst = 1.5;
Bconst = 2;
fxyt = @(tsol1, xsol1, ysol1) Aconst - (Bconst * xsol1) + xsol1^2 * ysol1 - xsol1;
gxyt = @(tsol1, xsol1, ysol1) Bconst * xsol1 - xsol1^2 * ysol1;
for i = 1:(length(tsol1)-1)
    k0 = fxyt(tsol1(i), xsol1(i), ysol1(i));
    l0 = gxyt(tsol1(i), xsol1(i), ysol1(i));
    k1 = fxyt(tsol1(i)+0.5*step,xsol1(i)+0.5*step*k0,ysol1(i)+0.5*step*l0);
    l1 = gxyt(tsol1(i)+0.5*step,xsol1(i)+0.5*step*k0,ysol1(i)+0.5*step*l0);
    k2 = fxyt((tsol1(i)+0.5*step),(xsol1(i)+0.5*step*k1),(ysol1(i)+0.5*step*l1));
    l2 = gxyt((tsol1(i)+0.5*step),(xsol1(i)+0.5*step*k1),(ysol1(i)+0.5*step*l1));
    k3 = fxyt((tsol1(i)+step),(xsol1(i)+k2*step),(ysol1(i)+l2*step));
    l3 = gxyt((tsol1(i)+step),(xsol1(i)+k2*step),(ysol1(i)+l2*step));
    xsol1(i+1) = xsol1(i) + (1/6)*step*(k0 + 2*k1 + 2*k2 + k3);
    ysol1(i+1) = ysol1(i) + (1/6)*step*(l0 + 2*l1 + 2*l2 + l3);
end
step = 0.1;
tsol2 = 0:step:tmax;
ysol2 = zeros(1, length(tsol2));
xsol2 = zeros(1, length(tsol2));
xsol2(1) = 1;
ysol2(1) = 1;
Aconst = 1.5;
Bconst = 2;
fxyt = @(tsol2, xsol2, ysol2) Aconst - (Bconst * xsol2) + xsol2^2 * ysol2 - xsol2;
gxyt = @(tsol2, xsol2, ysol2) Bconst * xsol2 - xsol2^2 * ysol2;
for i = 1:(length(tsol2)-1)
    k0 = fxyt(tsol2(i), xsol2(i), ysol2(i));
    l0 = gxyt(tsol2(i), xsol2(i), ysol2(i));
    k1 = fxyt(tsol2(i)+0.5*step,xsol2(i)+0.5*step*k0,ysol2(i)+0.5*step*l0);
    l1 = gxyt(tsol2(i)+0.5*step,xsol2(i)+0.5*step*k0,ysol2(i)+0.5*step*l0);
    k2 = fxyt((tsol2(i)+0.5*step),(xsol2(i)+0.5*step*k1),(ysol2(i)+0.5*step*l1));
    l2 = gxyt((tsol2(i)+0.5*step),(xsol2(i)+0.5*step*k1),(ysol2(i)+0.5*step*l1));
    k3 = fxyt((tsol2(i)+step),(xsol2(i)+k2*step),(ysol2(i)+l2*step));
    l3 = gxyt((tsol2(i)+step),(xsol2(i)+k2*step),(ysol2(i)+l2*step));
    xsol2(i+1) = xsol2(i) + (1/6)*step*(k0 + 2*k1 + 2*k2 + k3);
    ysol2(i+1) = ysol2(i) + (1/6)*step*(l0 + 2*l1 + 2*l2 + l3);
end
step = 1;
tsol3 = 0:step:tmax;
ysol3 = zeros(1, length(tsol3));
xsol3 = zeros(1, length(tsol3));
xsol3(1) = 1;
ysol3(1) = 1;
Aconst = 1.5;
Bconst = 2;
fxyt = @(tsol3, xsol3, ysol3) Aconst - (Bconst * xsol3) + xsol3^2 * ysol3 - xsol3;
gxyt = @(tsol3, xsol3, ysol3) Bconst * xsol3 - xsol3^2 * ysol3;
for i = 1:(length(tsol3)-1)
    k0 = fxyt(tsol3(i), xsol3(i), ysol3(i));
    l0 = gxyt(tsol3(i), xsol3(i), ysol3(i));
    k1 = fxyt(tsol3(i)+0.5*step,xsol3(i)+0.5*step*k0,ysol3(i)+0.5*step*l0);
    l1 = gxyt(tsol3(i)+0.5*step,xsol3(i)+0.5*step*k0,ysol3(i)+0.5*step*l0);
    k2 = fxyt((tsol3(i)+0.5*step),(xsol3(i)+0.5*step*k1),(ysol3(i)+0.5*step*l1));
    l2 = gxyt((tsol3(i)+0.5*step),(xsol3(i)+0.5*step*k1),(ysol3(i)+0.5*step*l1));
    k3 = fxyt((tsol3(i)+step),(xsol3(i)+k2*step),(ysol3(i)+l2*step));
    l3 = gxyt((tsol3(i)+step),(xsol3(i)+k2*step),(ysol3(i)+l2*step));
    xsol3(i+1) = xsol3(i) + (1/6)*step*(k0 + 2*k1 + 2*k2 + k3);
    ysol3(i+1) = ysol3(i) + (1/6)*step*(l0 + 2*l1 + 2*l2 + l3);
end

% calculating absolute average log oferror between different steps
errorx = log10(mean(abs(xsol(1:10:end)-xsol0)));
errorx0 = log10(mean(abs(xsol(1:100:end)-xsol1)));
errorx1 = log10(mean(abs(xsol(1:1000:end)-xsol2)));
errorx2 = log10(mean(abs(xsol(1:10000:end)-xsol3)));

% calculating best fit
coefficients = polyfit([1,2,3,4], [errorx, errorx0, errorx1, errorx2],1);
a = coefficients(1); %gradient estimate using least-squares

% plotting error and best fit
hold on;
plot([1;2;3;4], [errorx;errorx0;errorx1;errorx2], 'o', 'DisplayName', 'Average log error');
plot(1:0.1:4, polyval(coefficients,1:0.1:4), '-', 'DisplayName', 'Linear best fit');
title('RK4 difference in error as step size increases','FontSize',14);
xlabel('Step increased by h_c * 10^x', 'FontSize',14);
ylabel('log_{10} averaged error difference with h_c', 'FontSize',14);
grid on;
legend('location', 'northwest', 'FontSize',20);
set(gca,'FontSize',20);
hold off;