% Nathan Flynn
% 12/04/2022
% ME 573 - Final Project part 2
clc; close all; clear
format long;

%% Simulation inputs 
dX = 0.05;
dY = 0.05; 
L = 1;
gamma = 0;
nu = 0.01;
rho = 1;
C = 1;
tFinal  = 4;

%% Calcluated inputs
dT = C * min(1/4 * dX^2/nu);
time = 0:dT:tFinal;

%% Initializing nodes and applying initial conditions
[u,v,xu,yu,xv,yv,x_p,y_p,I,J] = generateNodes(dX, dY, L);
g = zeros(length(2:I),length(2:J));

% Main loop
for t = 1:length(0:dT:tFinal)
    %% Applying boundary conditions to the velocity field
    % Left boundary
    u(1,:) = 0;
    v(1,:) = -v(2,:);
    % Right boundary
    u(I,:) = 0;
    v(I+1,:) = -v(I,:);
    % Bottom boundary
    u(:,1) = -u(:,2);
    v(:,1) = 0;
    % Top boundary
    u(:,J+1) = 0;
    U = (u(:,J) + u(:,J+1))/2;
    u(:,J+1) = 2*U-u(:,J);
    v(:,J) = 0;

    %% Solving for velocity u and v
    [u, v] = solveUV(u,v,dX,dY,dT,I,J,gamma,nu);

    %% Computing the PPE source term
    for j = 2:J
        for i = 2:I
            g(i-1,j-1) = (rho/dT)*(((u(i,j) - u(i-1,j)))/dX + ((v(i,j) - v(i,j-1))/dY));
        end
    end
    nX = length(2:I);
    nY = length(2:J);
    B = reshape(g, [nX*nY, 1]);
    L = norm(B, 'inf');
end
%% Plots-------------------------------------------------------------------
% Plot 1 - u velocity
figure("units","normalized","position",[0,0.33,0.3,0.3])
surf(x_u,y_u,u)
xlabel('x')
ylabel('y')
set(gca('fontsize'),26)
title('u')
% Plot 2 - v velocity
figure('units','normalized','position',[0,0.01,0.3,0.3])
surf(x_v,y_v,v)
xlabel('x')
ylabel('y')
set(gca,'fontsize',26)
title('v')
% Plot 3 - pressure
figure('uints','normalized','position',[0.33,0.01,0.32,0.32])
surf(x_p,y_p,p)
set(gca,'fontsize',26)
xlabel('x')
ylabel('y')
title('Pressure')
% Plot 4 - validation (add this in later)
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
[nXu, nYu] = size(u); % Grabs the size of the rows and columns in u
[nXv, nYv] = size(v); % Grabs the size of the rows and columns in v
x_p = dX/2:dX:L-dX/2;
y_p = dY/2:dY:L-dY/2;
I = max([nXu nXv]);
J = max([nYu, nYv]);
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
u(2:I-1, 2:J) = u(2:I-1, 2:J) + dT*(nu * part_a_u - part_b_u - part_c_u);
v(2:I, 2:J-1) = v(2:I, 2:J-1) + dT*(nu * part_a_v - part_b_v - part_c_v);
end
%% Approximate functions for u
function out = d2u_dx2_plus_d2u_dy2(u,dX,dY,I,J)
nX = length(2:I-1);
nY = length(2:J);
d2u_dx2 = zeros(nX,nY);
d2u_dy2 = zeros(nX,nY);
for j = 2:J
    for i = 2:I-1
        d2u_dx2(i-1,j-1) = (u(i+1,j) - 2*u(i,j) + u(i-1,j))/dX^2;
        d2u_dy2(i-1,j-1) = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dY^2;
    end
end
out = d2u_dx2 + d2u_dy2;
end
function out = du2_dx(u,dX,I,J,gamma)
nX = length(2:I-1);
nY = length(2:J);
out = zeros(nX,nY);
for j = 2:J
    for i = 2:I-1
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
function out = duv_dy(u,v,dY,I,J,gamma)
nX = length(2:I-1);
nY = length(2:J);
out = zeros(nX,nY);
for j = 2:J
    for i = 2:I-1
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
function out = d2v_dx2_plus_d2v_dy2(v,dX,dY,I,J)
nX = length(2:I);
nY = length(2:J-1);
d2v_dx2 = zeros(nX,nY);
d2v_dy2 = zeros(nX,nY);
for j = 2:J-1
    for i = 2:I
        d2v_dx2(i-1,j-1) = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/dX^2;
        d2v_dy2(i-1,j-1) = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dY^2;
    end
end
out = d2v_dx2 + d2v_dy2;
end
function out = dv2_dy(v,dY,I,J,gamma)
nX = length(2:I);
nY = length(2:J-1);
out = zeros(nX,nY);
for j = 2:J-1
    for i = 2:I
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
function out = duv_dx(u,v,dX,I,J,gamma)
nX = length(2:I);
nY = length(2:J-1);
out = zeros(nX,nY);
for j = 2:J-1
    for i = 2:I
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