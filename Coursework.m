clear all
tmax = 1;
step = 0.5;
nsteps = round(tmax/step);
ysol = zeros(1, nsteps);
tsol = zeros(1, nsteps);
xsol = zeros(1, nsteps);
ysol(1) = 1;
xsol(1) = 0;
tsol(1) = 0;
Aconst = 2;
Bconst = 6;
for i = 2:nsteps+1
    k0 = step * fxy(tsol(i-1), xsol(i-1), ysol(i-1));
    l0 = step * gxy(tsol(i-1), xsol(i-1), ysol(i-1));
    k1 = step * fxy((tsol(i-1)+0.5*step), (xsol(i-1)+0.5*k0), (ysol(i-1)+0.5*l0));
    l1 = step * gxy((tsol(i-1)+0.5*step), (xsol(i-1)+0.5*k0), (ysol(i-1)+0.5*l0));
    k2 = step * fxy((tsol(i-1)+0.5*step), (xsol(i-1)+0.5*k1), (ysol(i-1)+0.5*l1));
    l2 = step * gxy((tsol(i-1)+0.5*step), (xsol(i-1)+0.5*k1), (ysol(i-1)+0.5*l1));
    k3 = step * fxy((tsol(i-1)+step), (xsol(i-1)+k2), (ysol(i-1)+l2));
    l3 = step * gxy((tsol(i-1)+step), (xsol(i-1)+k2), (ysol(i-1)+l2));
    xsol(i) = ysol(i-1) + (1/6)*(k0 + 2*k1 + 2*k2 + k3);
    ysol(i) = ysol(i-1) + (1/6)*(l0 + 2*l1 + 2*l2 + l3);
    tsol(i) = tsol(i-1) + step;

end
h=0.5; % step size 
x = 0:h:1; % Calculates upto y(1) 
y = zeros(1,length(x)); 
z = zeros(1,length(x)); 
y(1) = 0; % initial condition 
z(1) = 1; % initial condition
F_xyz = @(x,y,z) 2 -6*y + y^2 * z - y; % change the function as you desire 
G_xyz = @(x,y,z) 6 * y - y^2 * z; 
for i=1:(length(x)-1) % calculation loop 
    k_1 = F_xyz(x(i),y(i),z(i)); 
    L_1 = G_xyz(x(i),y(i),z(i)); 
    k_2 = F_xyz(x(i)+0.5*h,y(i)+0.5*h*k_1,z(i)+0.5*h*L_1); 
    L_2 = G_xyz(x(i)+0.5*h,y(i)+0.5*h*k_1,z(i)+0.5*h*L_1); 
    k_3 = F_xyz((x(i)+0.5*h),(y(i)+0.5*h*k_2),(z(i)+0.5*h*L_2)); 
    L_3 = G_xyz((x(i)+0.5*h),(y(i)+0.5*h*k_2),(z(i)+0.5*h*L_2)); 
    k_4 = F_xyz((x(i)+h),(y(i)+k_3*h),(z(i)+L_3*h)); % Corrected 
    L_4 = G_xyz((x(i)+h),(y(i)+k_3*h),(z(i)+L_3*h)); 
    y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; % main equation 
    z(i+1) = z(i) + (1/6)*(L_1+2*L_2+2*L_3+L_4)*h; % main equation 
end 
hold on
plot(x, z,'DisplayName', 'another one');
plot(tsol, ysol, 'DisplayName', 'Numerical sol')
legend('show')


function [fsol] = fxy(tval, xval, yval)
Aconst = 2;
Bconst = 6;
fsol = Aconst - Bconst * xval + (xval^2) * yval - xval;
end
function [gsol] = gxy(tval, xval, yval)
Bconst = 6;
gsol = Bconst * xval - (xval^2) * yval;
end