% 2D XY model, using time evolution method to calculate physical observable
% of finite tempature 

clear;
% close all;
clc;
format long
tic;

myseed = 2;
rng(myseed)

T_max = 1000;
dt = 1e-2;
M = 100;
t = 0:M*dt:T_max;
nt = length(t);

L = 100;
x = 1:L;
nx = length(x);
J = -3;
eta = 0.1;
theta = 1;

phi = zeros(L,L,2);
% phi(:,:,1) = 2*pi*ones(nx);
phi(:,:,1) = 2*pi*rand(nx);
vortex = zeros(L,L,nt);

titlename = strcat('\eta = ',num2str(eta), ',\theta = ',num2str(theta));
figure('Name',titlename);
set(gcf, 'position', [250 70 1600 800]);
len = 1;
u = cos(phi(:,:,1));
v = sin(phi(:,:,1));
xx = x.*ones(1,L) - u/2;
yy = x'.*ones(L,1) - v/2;
subplot(1,2,1)
h1 = quiver(xx,yy,u,v,0.6,'color','r','MaxHeadSize',0.7);
subplot(1,2,2)
vortex_p = zeros(L+1);
vortex_p(1:L,1:L) = vortex(:,:,1);
h2 = pcolor(vortex_p);

count = 2;
for i = 2:round(T_max/dt)+1
    phi = Heun_step(phi, dt, J, eta, theta);
%     phi = myFTCS(phi, dt, J, eta, 0);
    if mod(i+1,M) == 0
        vortex(:,:,count) = cal_vor(phi(:,:,1));
        u = cos(phi(:,:,1));
        v = sin(phi(:,:,1));
        xx = x.*ones(1,L) - u/2;
        yy = x'.*ones(L,1) - v/2;
        h1.XData = xx;
        h1.YData = yy;
        h1.UData = u;
        h1.VData = v;
        vortex_p(1:L,1:L) = vortex(:,:,count);
        h2.CData = vortex_p;
        drawnow
        pause(0.1)
        count = count + 1;
    end
end

toc;

function y = myFTCS(phi, dt, J, eta, field)
c1 = f(phi, J, eta, field);
c2 = f(phi+c1*dt/2, J, eta, field);
c3 = f(phi+c2*dt/2, J, eta, field);
c4 = f(phi+c3*dt, J, eta, field);
y = phi + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = Heun_step(phi, dt, J, eta, theta)
L = length(phi);
field = randn(L)*sqrt(2*eta*theta/dt);
z = phi + dt*f(phi, J, eta, field);
y = phi + dt*f((phi+z)/2, J, eta, field);
end

function y = f(x, J, eta, field)
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
y(:,:,2) = -eta*dphi + J*tau + field;
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

% function y = cal_vor(phi)
% diff_l = phi-circshift(phi,1,1);
% diff_l = atan2(sin(diff_l),cos(diff_l));
% diff_r = -circshift(diff_l,-2,1);
% diff_u = phi-circshift(phi,1,2);
% diff_u = atan2(sin(diff_u),cos(diff_u));
% diff_d = -circshift(diff_u,-2,2);
% y = diff_l + diff_r + diff_u + diff_d;
% y = y/pi;
% end

% function y = cal_vor(phi)
% phi = phi/pi;
% diff_l = circshift(phi,1,1);
% diff_r = circshift(phi,-1,1);
% diff_u = circshift(phi,1,2);
% diff_d = circshift(phi,-1,2);
% d1 = mod(diff_l - diff_u + 1,2)-1;
% d2 = mod(diff_u - diff_r + 1,2)-1;
% d3 = mod(diff_r - diff_d + 1,2)-1;
% d4 = mod(diff_d - diff_l + 1,2)-1;
% y = d1+d2+d3+d4;
% end

function y = cal_vor(phi)
phi = phi/pi;
diff1 = circshift(phi,1,1);
diff2 = circshift(phi,[1,1]);
diff3 = circshift(phi,1,2);
d1 = mod(phi - diff1 + 1,2)-1;
d2 = mod(diff1 - diff2 + 1,2)-1;
d3 = mod(diff2 - diff3 + 1,2)-1;
d4 = mod(diff3 - phi + 1,2)-1;
y = d1+d2+d3+d4;
end