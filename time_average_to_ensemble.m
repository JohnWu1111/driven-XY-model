% 2D XY model, using time evolution method to calculate physical observable
% of finite tempature 

clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

T_pre = 1000;
dt = 1e-2;
t_pre = 0:dt:T_pre;
nt_pre = length(t_pre);
n_mea = 1000;
M = 100;
L = 50;
x = 1:L;
nx = length(x);
J = -1;
eta = 0.1;
theta = 0.25;
mass = 1;

phi = zeros(L,L,2);
phi(:,:,1) = 2*pi*rand(nx);
order = zeros(1,nt_pre);
order(1) = sum(cos(phi(:,:,1)),'all');
corre = zeros(L,1);

% pre-thermalization
for i = 2:nt_pre
    phi = Heun_step(phi, dt, J, eta, theta, mass);
%     phi = myFTCS(phi, dt, J, eta, 0);
    order(i) = sum(cos(phi(:,:,1)),'all');
end
order = order/L^2;

% measurement
for i = 1:M*n_mea
    phi = Heun_step(phi, dt, J, eta, theta, mass);
%     phi = myFTCS(phi, dt, J, eta, 0);
    if mod(i,M) == 0
        corre = corre + cal_corre(phi(:,:,1));
    end
end
corre = corre/n_mea;

corre_p = circshift(corre,1);
factor = -theta/(2*pi*J);

toc;

titlename = strcat('L = ',num2str(L),',J = ',num2str(J),',\eta = ',num2str(eta), ',\theta = ',num2str(theta));
figure('Name',titlename);
set(gcf, 'position', [250 70 1600 900]);
subplot(1,3,1)
plot(t_pre,order)
subplot(1,3,2)
loglog(x(1:L/2),abs(corre_p(1:L/2)))
subplot(1,3,3)
semilogy(x(1:L/2),abs(corre_p(1:L/2)))

function y = myFTCS(phi, dt, J, eta, field) %#ok<DEFNU>
c1 = f(phi, J, eta, field);
c2 = f(phi+c1*dt/2, J, eta, field);
c3 = f(phi+c2*dt/2, J, eta, field);
c4 = f(phi+c3*dt, J, eta, field);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = Heun_step(phi, dt, J, eta, theta, mass)
L = length(phi);
field = randn(L)*sqrt(2*eta*theta/dt);
z = phi + dt*f(phi, J, eta, field, mass);
y = phi + dt*f((phi+z)/2, J, eta, field, mass);
end

function y = f(x, J, eta, field, mass)
y = x;
phi = x(:,:,1);
dphi = x(:,:,2);
y(:,:,1) = dphi;

diff_l = sin(phi-circshift(phi,1,1));
diff_u = sin(phi-circshift(phi,1,2));
diff_r = -circshift(diff_l,-1,1);
diff_d = -circshift(diff_u,-1,2);
tau = diff_l + diff_r + diff_u + diff_d;
% tau = sin(phi-circshift(phi,1,1)) + sin(phi-circshift(phi,-1,1)) + sin(phi-circshift(phi,1,2)) + sin(phi-circshift(phi,-1,2));
y(:,:,2) = (-eta*dphi + J*tau + field)/mass;
end

function y = cal_corre(phi)
L = length(phi);
y = zeros(L,1);
for i = 1:L
    diff_l = cos(phi-circshift(phi,i,1));
    diff_r = circshift(diff_l,-i,1);
    diff_u = cos(phi-circshift(phi,i,2));
    diff_d = circshift(diff_u,-i,2);
    y(i) = sum(diff_l + diff_r + diff_u + diff_d,'all')/(4*L^2);
%     y(i) = sum(cos(phi-circshift(phi,i,1)) + cos(phi-circshift(phi,-i,1)) + cos(phi-circshift(phi,i,2)) + cos(phi-circshift(phi,-i,2)),'all')/(4*L^2);
end
end