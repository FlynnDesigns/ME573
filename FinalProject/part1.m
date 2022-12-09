% Nathan Flynn
% 12/04/2022
% ME 573 - Final Project
clc; close all; clear
format long;
%(i)
nu = 0.5;
[x_p,y_p,uCenter,vCenter,umag,time_i,KE_i] = run_sim(nu);
% Plot 1
figure('units','normalized','position',[0.3,0.65,0.4,0.4])
quiver(x_p,y_p,uCenter,vCenter,2)
set(gca,'fontsize',26)
xlabel('x_p')
ylabel('y_p')
% xlim([-0.125 1.125])
% ylim([-0.125 1.125])
[~,b] = title('Velocity Field \nu=0.5','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
% Plot 2
figure('units','normalized','position',[0.3,0,0.4,0.4])
contour(x_p,y_p,umag,'ShowText','on')
set(gca,'fontsize',26)
xlabel('x_p')
ylabel('y_p')
% xlim([-0.125 1.125])
% ylim([-0.125 1.125])
[~,b] = title('Velocity Magnitude \nu=0.5','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
%(ii)
nu = 0.05;
[x_p,y_p,uCenter,vCenter,umag,time_ii,KE_ii] = run_sim(nu);
% Plot 1
figure('units','normalized','position',[0.3,0.65,0.4,0.4])
quiver(x_p,y_p,uCenter,vCenter,2)
set(gca,'fontsize',26)
xlabel('x_p')
ylabel('y_p')
% xlim([-0.125 1.125])
% ylim([-0.125 1.125])
[~,b] = title('Velocity Field \nu=0.05','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
% Plot 2
figure('units','normalized','position',[0.3,0,0.4,0.4])
contour(x_p,y_p,umag,'ShowText','on')
set(gca,'fontsize',26)
xlabel('x_p')
ylabel('y_p')
% xlim([-0.125 1.125])
% ylim([-0.125 1.125])
[~,b] = title('Velocity Magnitude \nu=0.05','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
%(iii)
nu = 0.005;
[x_p,y_p,uCenter,vCenter,umag,time_iii,KE_iii] = run_sim(nu);
% Plot 1
figure('units','normalized','position',[0.3,0.65,0.4,0.4])
quiver(x_p,y_p,uCenter,vCenter,2)
set(gca,'fontsize',26)
xlabel('x_p')
ylabel('y_p')
% xlim([-0.125 1.125])
% ylim([-0.125 1.125])
[~,b] = title('Velocity Field \nu=0.005','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
% Plot 2
figure('units','normalized','position',[0.3,0,0.4,0.4])
contour(x_p,y_p,umag,'ShowText','on')
set(gca,'fontsize',26)
xlabel('x_p')
ylabel('y_p')
% xlim([-0.125 1.125])
% ylim([-0.125 1.125])
[~,b] = title('Velocity Magnitude \nu=0.005','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
% Plot 3
figure('units','normalized','position',[0.65,0.01,0.5,0.5])
hold on;
plot(time_i,KE_i,'-')
plot(time_ii,KE_ii)
plot(time_iii,KE_iii)
legend('\nu = 0.5','\nu=0.05','\nu=0.005')
set(gca,'fontsize',26)
xlabel('time (s)')
[~,b] = title('Total Kinetic Energy','\Deltax=0.05, \Deltay=0.05, \gamma=0.5, C=0.005, Tfinal=1');
b.FontSize = 16;
%% Function to run the simulation (used to generate multiple plots)
function [x_p,y_p,uCenter,vCenter,umag,time,KE] = run_sim(nu)
% Simulation inputs 
dX = 0.05;
dY = 0.05; 
L = 1;
gamma = 0.5;
C = 0.005;
tFinal  = 1;

% Calcluated inputs
dT = C * min(1/4 * dX^2/nu);
time = 0:dT:tFinal;

% Initializing nodes and applying initial conditions
[u,v,~,~,~,~,x_p,y_p,I,J] = generateNodes(dX, dY, L);

% Allocating arrays
KE = zeros(1, length(time));
uCenter = zeros(I,J);
vCenter = zeros(I,J);

% Solving u and v
for i = 1:length(0:dT:tFinal)
    %% Applying boundary conditions to the velocity field
    % Left boundary
    u(1,2:J) = 0;
    % Right boundary
    u(I,2:J) = 0;
    % Bottom boundary
    u(2:I,1) = 0;
    % Top boundary
    u(2:I,J+1) = 0;
    %% Solving for u and v
    [u, v] = solveUV(u,v,dX,dY,dT,I,J,gamma,nu);

    %% Computing the umag and kenetic energy
    % Getting the centers of all of the u nodes
    for k = 2:I
        for  j = 2:J+1
            uCenter(k,j-1) = (u(k,j) + u(k,j-1))/2;
        end
    end
    for k = 2:J+1
        for j = 2:I
            vCenter(k-1,j) = (v(k,j) + v(k-1,j))/2;
        end
    end
    % Calculating the velocity magnitude
    umag = sqrt(uCenter.^2 + vCenter.^2);
     % Calculating the kenetic energy (KE)
    KE(i) = KE(i) + 0.5*dX*dY*sum((uCenter.^2 + vCenter.^2),"All");
end
end
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
u = u';
v = v';

% Points for pressure 
x_p = 0:dX:L;
y_p = 0:dY:L;

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
        d2u_dx2(i,j) = (u(i+1,j) - 2*u(i,j) + u(i-1,j))/dX^2; %Checked
        d2u_dy2(i,j) = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dY^2; %Checked
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
            *((abs(u(i,j)+u(i+1,j))/2) ...
            *((u(i,j)-u(i+1,j))/2) ...
            -(abs(u(i-1,j)+u(i,j))/2) ...
            *((u(i-1,j)-u(i,j))/2)); %Checked
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
            *(abs((v(i,j) + v(i+1,j))/2) ...
            *((u(i,j) - u(i,j+1))/2) ...
            -(abs((v(i,j-1) + v(i+1,j-1))/2)) ...
            *((u(i,j-1) - u(i,j))/2)); %Checked
     end
end
end
%% Approximate functions for v
function out = d2v_dx2_plus_d2v_dy2(v,dX,dY,I,J)
d2v_dx2 = zeros(I+1,J);
d2v_dy2 = zeros(I+1,J);
for j = 2:J-1
    for i = 2:I
        d2v_dx2(i,j) = (v(i+1,j) - 2*v(i,j) + v(i-1,j))/dX^2; %Checked
        d2v_dy2(i,j) = (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dY^2; %Checked
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
            -((v(i,j-1) + v(i,j))/2)^2) ...
            +(gamma/dY)...
            *((abs(v(i,j)+v(i,j+1))/2) ...
            *((v(i,j)-v(i,j+1))/2) ...
            -(abs(v(i,j-1)+v(i,j))/2) ...
            *((v(i,j-1)-v(i,j))/2)); %Checked
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
            *((abs((u(i,j+1)+u(i,j))/2)) ...
            *((v(i,j) - v(i+1,j))/2) ...
            -(abs((u(i-1,j+1) + u(i-1,j))/2)) ...
            *((v(i-1,j) - v(i,j))/2)); %Checked
     end
end
end