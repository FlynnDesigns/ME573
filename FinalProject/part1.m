% Nathan Flynn
% 12/04/2022
% ME 573 - Final Project
clc; close all; clear
format long;

% Simulation inputs 
dX = 0.05;
dY = 0.05; 
L = 1;
gamma = 0.5;
nu = 0.5;
C = 0.01;
tFinal  = 1.2;

% Calcluated inputs
dT = C * min(1/4 * dX^2/nu);
time = 0:dT:tFinal;

% Initializing nodes and applying initial conditions
[u,v,xu,yu,xv,yv,x_p,y_p,I,J] = generateNodes(dX, dY, L);

% Allocating arrays
KE = zeros(1, length(time));
uCenter = zeros(I-1,J-1);
vCenter = zeros(I-1,J-1);

% Solving u and v
for i = 1:length(0:dT:tFinal)
    %% Applying boundary conditions to the velocity field
    % Left boundary
    u(1,2:J) = 0;
    v(1,2:J) = -v(2,2:J);
    % Right boundary
    u(I,:) = 0;
    v(I+1,2:J) = -v(I,2:J);
    % Bottom boundary
    u(:,1) = 0;
    v(2:I,1) = 0;
    % Top boundary
    u(:,J+1) = 0;
    v(2:I,J) = 0;
   
    %% Solving for u and v
    [u, v] = solveUV(u,v,dX,dY,dT,I,J,gamma,nu);

    %% Computing the umag and kenetic energy
    % Getting the centers of all of the u nodes
    for j = 2:J
        for  k = 2:I
            uCenter(k-1,j-1) = (u(k-1,j) + u(k,j))/2;
            vCenter(k-1,j-1) = (v(k,j-1) + v(k,j))/2;
        end
    end
    % Calculating the velocity magnitude
    umag = sqrt(uCenter.^2 + vCenter.^2);
    % Calculating the kenetic energy (KE)
    KE(i) = KE(i) + sum(0.5*dX*dY*(uCenter.^2 + vCenter.^2),"All");
end
%% Plots-------------------------------------------------------------------
% Plot 1
figure('units','normalized','position',[0.3,0.65,0.4,0.4])
quiver(x_p,y_p,uCenter,vCenter,2)
set(gca,'fontsize',26)
xlim([0 1])
ylim([0 1])
xlabel('x_p')
ylabel('y_p')
title('Velocity Field')
% Plot 2
figure('units','normalized','position',[0.3,0,0.4,0.4])
contour(x_p,y_p,umag,'ShowText','on')
set(gca,'fontsize',26)
xlim([0 1])
ylim([0 1])
xlabel('x_p')
ylabel('y_p')
title('Velocity Magnitude')
% Plot 3
figure('units','normalized','position',[0.65,0.01,0.5,0.5])
plot(time,KE,'-')
set(gca,'fontsize',26)
xlabel('time (s)')
title('Total Kinetic Energy')
%% Function to generate U nodes
function [u,v,xu,yu,xv,yv,x_p,y_p,I,J] = generateNodes(dX, dY, L)
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

% Points for pressure 
x_p = dX/2:dX:L-dX/2;
y_p = dY/2:dY:L-dY/2;
% Indicies 
I = length(xu);
J = length(yv);
end
%% Function to solve for u and v
function [u, v] = solveUV(u,v,dX,dY,dT,I,J,gamma,nu)
% Part a
part_a_u = d2u_dx2_plus_d2u_dy2(u,dX,dY,I,J);
part_a_v = d2v_dx2_plus_d2v_dy2(v,dX,dY,I,J);
% Part b
part_b_u = du2_dx(u,dX,I,J,gamma);
part_b_v = dv2_dy(v,dY,I,J,gamma);
% Part c
part_c_u = duv_dy(u,v,dY,I,J,gamma);
part_c_v = duv_dx(u,v,dX,I,J,gamma);
% Computing u and v
u = u + dT*(nu * part_a_u - part_b_u - part_c_u);
v = v + dT*(nu * part_a_v - part_b_v - part_c_v);
end
%% Approximate functions for u
function out = d2u_dx2_plus_d2u_dy2(u,dX,dY,I,J)
d2u_dx2 = zeros(I,J+1);
d2u_dy2 = zeros(I,J+1);
for j = 2:J
    for i = 2:I-1
        d2u_dx2(i,j) = (u(i+1,j) - 2*u(i,j) + u(i-1,j))/dX^2;
        d2u_dy2(i,j) = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dY^2;
    end
end
out = d2u_dx2 + d2u_dy2;
end
function out = du2_dx(u,dX,I,J,gamma)
out = zeros(I,J+1);
for j = 2:J
    for i = 2:I-1
        out(i,j) = (1/dX)...
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
function out = duv_dy(u,v,dY,I,J,gamma)
out = zeros(I,J+1);
for j = 2:J
    for i = 2:I-1
        out(i,j) = (1/dY)...
            *(((v(i,j)+v(i+1,j))/2)...
            *((u(i,j+1)+u(i,j))/2)...
            -((v(i,j-1)+v(i+1,j-1))/2)...
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
function out = d2v_dx2_plus_d2v_dy2(v,dX,dY,I,J)
d2v_dx2 = zeros(I+1,J);
d2v_dy2 = zeros(I+1,J);
for j = 2:J-1
    for i = 2:I
        d2v_dx2(i,j) = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/dX^2;
        d2v_dy2(i,j) = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dY^2;
    end
end
out = d2v_dx2 + d2v_dy2;
end
function out = dv2_dy(v,dY,I,J,gamma)
out = zeros(I+1,J);
for j = 2:J-1
    for i = 2:I
        out(i,j) = (1/dY)...
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
function out = duv_dx(u,v,dX,I,J,gamma)
out = zeros(I+1,J);
for j = 2:J-1
    for i = 2:I
        out(i,j) = (1/dX)...
            *(((v(i+1,j)+v(i,j))/2)...
            *((u(i,j+1)+u(i,j))/2)...
            -((v(i-1,j)+v(i,j))/2)...
            *((u(i-1,j+1)+u(i-1,j))/2))...
            +(gamma/dX)...
            *(abs((u(i,j+1)+u(i,j)))/2 ...
            *(v(i,j) - v(i+1,j))/2 ...
            -(abs((u(i-1,j+1) + u(i-1,j)))/2) ...
            *(v(i-1,j) - v(i,j))/2);
     end
end
end