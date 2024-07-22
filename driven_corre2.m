% 2D XY model, using time evolution method to calculate physical observable
% of finite tempature with staggered interaction strength as a drive

clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

T_max = 2000;
M = 100;
dt = 1e-2;
t = 0:M*dt:T_max;
nt = length(t);
pre = 1000;

L = 50;
x = 1:L;
nx = length(x);
J = -1;
eta = 0.1;
theta = 1;
mass = 1;

period = 1;
% driven = [2,3];
% driven = [0.5,2];
% driven = [1/3,2/3];
driven = [1/2,1/4];
% driven = [1,0.9];
% driven = [0.8,0.9];
% driven = [0.8,1.5];
% driven = [0.2,0.9];

phi = zeros(L,L,2);
% phi(:,:,1) = 2*pi*rand(nx);
order = zeros(1,nt);
order(1) = sum(cos(phi(:,:,1)),'all');
corre = zeros(L,1);
corre_mid = zeros(L,1);
E = zeros(1,nt);
E(1) = J*cal_energy(phi(:,:,1));

count = 2;
for i = 2:round(T_max/dt)+1
    temp = mod(floor(i*dt/period),2)+1;
    phi = Heun_step(phi, dt, J, eta, driven(temp),mass);
%     phi = myFTCS(phi, dt, J, eta, 0);
    if mod(i+1-M/2,M) == 0
        if count > pre
            corre_mid = corre_mid + cal_corre(phi(:,:,1));
        end
    end
    if mod(i+1,M) == 0
        order(:,count) = sum(cos(phi(:,:,1)),'all');
        E(count) = J*cal_energy(phi(:,:,1));
        if count > pre
            corre = corre + cal_corre(phi(:,:,1));
        end
        count = count + 1;
    end
end
order = order/L^2;
corre = corre/(nt-pre);
corre_mid = corre_mid/(nt-pre);

% measurement
% for i = 1:M*n_mea
%     phi = Heun_step(phi, dt, J, eta, theta);
% %     phi = myFTCS(phi, dt, J, eta, 0);
%     if mod(i,M) == 0
%         corre = corre + cal_corre(phi(:,:,1));
%     end
% end
% corre = corre/n_mea;
% 
corre_p = circshift(corre,1);
corre_mid_p = circshift(corre_mid,1);

toc;

titlename = strcat('L = ',num2str(L),',\eta = ',num2str(eta), ',dt = ',num2str(dt), ',T1 = ',num2str(driven(1)), ',T2 = ',num2str(driven(2)), ',period = ',num2str(period), ',mass = ',num2str(mass));
figure('Name',titlename);
set(gcf, 'position', [250 70 1600 900]);
subplot(2,2,1)
plot(t,order)
subplot(2,2,2)
plot(t(pre:end),E(pre:end))
subplot(2,2,3)
loglog(x(1:L/2),abs(corre_p(1:L/2)),x(1:L/2),abs(corre_mid_p(1:L/2)))
legend()
subplot(2,2,4)
semilogy(x(1:L/2),abs(corre_p(1:L/2)))

function y = myFTCS(phi, dt, J, eta, field)
c1 = f(phi, J, eta, field);
c2 = f(phi+c1*dt/2, J, eta, field);
c3 = f(phi+c2*dt/2, J, eta, field);
c4 = f(phi+c3*dt, J, eta, field);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = Heun_step(phi, dt, J, eta, theta,mass)
L = length(phi);
field = randn(L)*sqrt(2*eta*theta/dt);
z = phi + dt*f(phi, J, eta, field,mass);
y = phi + dt*f((phi+z)/2, J, eta, field,mass);
end

function y = f(x, J, eta, field,mass)
y = x;
phi = x(:,:,1);
dphi = x(:,:,2);
y(:,:,1) = dphi;

diff_l = sin(phi-circshift(phi,1,1));
diff_r = -circshift(diff_l,-1,1);
diff_u = sin(phi-circshift(phi,1,2));
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

function y = cal_energy(phi)
L = length(phi);
diff_l = cos(phi-circshift(phi,1,1));
diff_r = circshift(diff_l,-1,1);
diff_u = cos(phi-circshift(phi,1,2));
diff_d = circshift(diff_u,-1,2);
y = sum(diff_l + diff_r + diff_u + diff_d,'all')/(2*L^2);
end