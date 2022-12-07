% Nathan Flynn
% 12/04/2022
% ME 573 - Final Project
clc; close all; clear
format long;

% Simulation inputs 
dX = 0.05;
dY = 0.05; 
L = 1;
gamma = 1;
nu = 0.005;
C = 0.0001;
tFinal  = 0.02;

% Calcluated inputs
dT = C * min(1/4 * dX^2/nu);
time = 0:dT:tFinal;

% Initializing nodes and applying initial conditions
[u,v,xu,yu,xv,yv,x_p,y_p,nXu,nYu,nXv,nYv] = generateNodes(dX, dY, L);

% Allocating arrays
KE = zeros(1, length(time));
KE(1) = 0;
vCenter = zeros(nYv-1, nYv-1);
uCenter = zeros(nXu-1,nXu-1);

% Solving u and v
for i = 1:length(0:dT:tFinal)
    % Applying boundary conditions
    u(1,:) = 0;
    v(1,:) = 0;
    u(nXu,:) = 0;
    v(nXu+1,:) = 0;
    v(:,1) = 0;
    v(:,nYv) = 0;
   
    % Solving for u and v
    [u, v] = solveUV(u, v, dX, dY, dT, nXu, nYu, nXv, nYv, 1, 1);

    % Getting the centers of all of the u nodes
    uCenter = zeros(nXu-1,nXu-1);
    for j = 2:nYu-1
        for  k = 2:nXu
            uCenter(k-1,j-1) = (u(k-1,j) + u(k,j))/2;
        end
    end

    % Getting the centers of all of the v nodes 
    vCenter = zeros(nYv-1, nYv-1);
    for j = 2:nYv % running through y
        for k = 2:nXv-1 % running through x
            vCenter(k-1,j-1) = (v(k,j-1) + v(k,j))/2;
        end
    end

    % Calculating the velocity magnitude
    umag = sqrt(uCenter.^2 + vCenter.^2);

    % Calculating the kenetic energy (KE)
    KE(i) = KE(i) + sum(0.5*dX*dY*(uCenter.^2 + vCenter.^2),"All");
end

%% Plot 1
figure('units','normalized','position',[0.3,0.65,0.4,0.4])
quiver(x_p,y_p,uCenter,vCenter,2)
set(gca,'fontsize',26)
xlim([-0.2,1.2])
ylim([0,1])
xlabel('x_p')
ylabel('y_p')
title('Velocity Field')
%% Plot 2
figure('units','normalized','position',[0.3,0,0.4,0.4])
contour(x_p,y_p,umag,'ShowText','on')
set(gca,'fontsize',26)
xlim([0,1])
ylim([0,1])
xlabel('x_p')
ylabel('y_p')
title('Velocity Magnitude')
%% Plot 3
figure('units','normalized','position',[0.65,0.01,0.5,0.5])
plot(time,KE,'-')
set(gca,'fontsize',26)
xlabel('time (s)')
%legend('\nu=0.005','\nu=0.01','\nu=0.05')
title('Total Kinetic Energy')
%% Function to generate U nodes
function [u,v,xu,yu,xv,yv,x_p,y_p,nXu,nYu,nXv,nYv] = generateNodes(dX, dY, L)
% Nodes in x
xu = 0:dX:L;
xv = -dX/2:dX:L+dX/2;

% Nodes in y
yu = -dY/2:dY:L+dY/2;
yv = 0:dY:L;

% Creating points
[y_u, x_u] = meshgrid(xu, yu);
[y_v, x_v] = meshgrid(xv, yv);

% Applying initial conditions
u = ((sin(pi*x_u).^2).*sin(2*pi*y_u));
v = -((sin(pi*y_v)).^2).*sin(2*pi*x_v);
u =u';
v = v';
[nXu, nYu] = size(u); % Grabs the size of the rows and columns in u
[nXv, nYv] = size(v); % Grabs the size of the rows and columns in v
x_p = dX/2:dX:L-dX/2;
y_p = dY/2:dY:L-dY/2;
end
%% Function to solve for u and v
function [u, v] = solveUV(u,v,dX,dY,dT,nXu,nYu,nXv,nYv,gamma,nu)
% Part a
part_a_u = d2u_dx2_plus_d2u_dy2(u,dX,dY,nXu,nYu);
part_a_v = d2v_dx2_plus_d2v_dy2(v,dX,dY,nXv,nYv);
% Part b
part_b_u = du2_dx(u,dX,nXu,nYu,gamma);
part_b_v = dv2_dy(v,dY,nXv,nYv,gamma);
% Part c
part_c_u = duv_dy(u,v,dY,nXu,nYu,gamma);
part_c_v = duv_dx(u,v,dX,nXv,nYv,gamma);
% Computing u and v
u(2:nXu-1, 2:nYu-1) = u(2:nXu-1, 2:nYu-1) + dT*(nu * part_a_u - part_b_u - part_c_u);
v(2:nXv-1, 2:nYv-1) = v(2:nXv-1, 2:nYv-1) + dT*(nu * part_a_v - part_b_v - part_c_v);
end
%% Approximate functions for u
function out = d2u_dx2_plus_d2u_dy2(u,dX,dY,nXu,nYu)
d2u_dx2 = zeros(nXu-2,nYu-2);
d2u_dy2 = zeros(nXu-2,nYu-2);
for j = 2:nYu-1
    for i = 2:nXu-1
        d2u_dx2(i-1,j-1) = (u(i+1,j) - 2*u(i,j) + u(i-1,j))/dX^2;
        d2u_dy2(i-1,j-1) = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dY^2;
    end
end
out = d2u_dx2 + d2u_dy2;
end
function out = du2_dx(u,dX,nXu,nYu,gamma)
out = zeros(nXu-2,nYu-2);
for j = 2:nYu-1
    for i = 2:nXu-1
        out(i-1,j-1) = (1/dX)...
            *(((u(i,j) + u(i+1,j))/2)^2 ...
            -(((u(i-1,j) + u(i,j))/2)^2)) ...
            +(gamma/dX)...
            *(abs(u(i,j)+u(i+1,j))/2 ...
            *(u(i,j)-u(i+1,j))/2 ...
            -abs(u(i-1,j)+u(i,j))/2 ...
            *(u(i-1,j)-u(i,j))/2);
    end
end
end
function out = duv_dy(u,v,dY,nXu,nYu,gamma)
out = zeros(nXu-2,nYu-2);
for j = 2:nYu-1
    for i = 2:nXu-1
        out(i-1,j-1) = (1/dY)...
            *(((v(i,j)+v(i+1,j))/2)...
            *((u(i,j+1)+u(i,j))/2)...
            -(v(i,j-1)+v(i+1,j-1)/2)...
            *((u(i,j-1)+u(i,j))/2))...
            +(gamma/dY)...
            *(abs((v(i,j) + v(i+1,j)))/2 ...
            *(u(i,j) - u(i,j+1))/2 ...
            -(abs((v(i,j-1) + v(i+1,j-1)))/2) ...
            *(u(i,j-1) - u(i,j))/2);
     end
end
end
%% Approximate functions for v
function out = d2v_dx2_plus_d2v_dy2(v,dX,dY,nXv,nYv)
d2v_dx2 = zeros(nXv-2,nYv-2);
d2v_dy2 = zeros(nXv-2,nYv-2);
for j = 2:nYv-1
    for i = 2:nXv-1
        d2v_dx2(i-1,j-1) = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/dX^2;
        d2v_dy2(i-1,j-1) = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dY^2;
    end
end
out = d2v_dx2 + d2v_dy2;
end
function out = dv2_dy(v,dY,nXv,nYv,gamma)
out = zeros(nXv-2,nYv-2);
for j = 2:nYv-1
    for i = 2:nXv-1
        out(i-1,j-1) = (1/dY)...
            *(((v(i,j) + v(i,j+1))/2)^2 ...
            -(((v(i,j-1) + v(i,j))/2)^2)) ...
            +(gamma/dY)...
            *(abs(v(i,j)+v(i,j+1))/2 ...
            *(v(i,j)-v(i,j+1))/2 ...
            -abs(v(i,j-1)+v(i,j))/2 ...
            *(v(i,j-1)-v(i,j))/2);
    end
end
end
function out = duv_dx(u,v,dX,nXv,nYv,gamma)
out = zeros(nXv-2,nYv-2);
for j = 2:nYv-1
    for i = 2:nXv-1
        out(i-1,j-1) = (1/dX)...
            *(((v(i+1,j)+v(i,j))/2)...
            *((u(i,j+1)+u(i,j))/2)...
            -(v(i-1,j)+v(i,j)/2)...
            *((u(i-1,j+1)+u(i-1,j))/2))...
            +(gamma/dX)...
            *(abs((u(i,j+1)+u(i,j)))/2 ...
            *(v(i,j) - v(i+1,j))/2 ...
            -(abs((u(i-1,j+1) + u(i-1,j)))/2) ...
            *(v(i-1,j) - v(i,j))/2);
     end
end
end