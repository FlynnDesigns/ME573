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
Tolfac = 10^(-7);
uLid = 1;

%% Calcluated inputs
dT = C * min(1/4 * dX^2/nu, dX/uLid);
time = 0:dT:tFinal;

%% Initializing nodes and applying initial conditions
[u,v,xu,yu,xv,yv,p_k,x_p,y_p,I,J] = generateNodes(dX, dY, L);
g = zeros(I+1,J+1);
res = zeros(I+1,J+1);
eW = zeros(1,I+1);
eE = zeros(1,I+1);
eN = zeros(1,J+1);
p_kp1 = p_k;
pois_iter = zeros(1, length(time));
pois_total = 0;

%% Calculating indicator functions
% eW
for i = 2:I
    if i == 2
        eW(i) = 0;
    end
    if i > 2
        eW(i) = 1;
    end
end
% eE
for i = 2:I
    if i < I
        eE(i) = 1;
    end
    if i == I 
        eE(i) = 0;
    end
end
% eN
for j = 2:J
    if j < J
        eN(j) = 1;
    end
    if j == J
        eN(j) = 0;
    end
end
% Main loop
for t = 1:length(0:dT:tFinal)
    %% Applying boundary conditions to the velocity field
    j = 2:J;
    i = 2:I;
    % Left boundary
    u(1,j) = 0; % Checked
    v(1,j) = -v(2,j); % Checked
    % Right boundary
    u(I,j) = 0; % Checked
    v(I+1,j) = -v(I,j); % Checked
    % Bottom boundary
    u(i,1) = -u(i,2); %Checked
    v(i,1) = 0; %Checked
    % Top boundary
    u(i,J+1) = 2*uLid - u(i,J); %Checked
    v(i,J) = 0; %Checked

    %% Applying boundary conditions to the pressure field
    % Left boundary
    i = 2:I;
    j = 2:J;
    p_k(1,j) = -p_k(2,j);
    % Right boundary
    p_k(end,j) = -p_k(end-1,j);
    % Bottom boundary
    p_k(i,1) = 0;
    % Top boundary
    p_k(i,end) = -p_k(i,end-1);

    %% Solving for velocity u and v
    [u_star, v_star] = solveUV(u,v,dX,dY,dT,I,J,gamma,nu);

    %% Computing the PPE source term
    for j = 2:J
        for i = 2:I
            g(i,j) = (rho/dT)*((u_star(i,j) - u_star(i-1,j))/dX + (v_star(i,j) - v_star(i,j-1))/dY);
        end
    end
    B = reshape(g, [(I+1)*(J+1), 1]);
    L = norm(B, 'inf');
    Tol = L*Tolfac;

    %% Calculating pressure
    L_inf_residual = Tol*10;
    iter = 0; 
    pois_total = 0; 
    while (L_inf_residual > Tol)
        for j = 2:J
            for i = 2:I
                % Caclulating spectral radius and optimum omega
                lambda = 0.5*(cos(pi*i)/I + cos(pi*j)/J);
                omega = 2/(1 + sqrt(1-lambda^2));

                % Caclulating p_kp1
                p_kp1(i,j) = p_k(i,j)*(1-omega) ...
                    + omega/(eE(i)+eW(i)+eN(j)+1) ...
                    * (((eE(i)*p_k(i+1,j)+eW(i)*p_kp1(i-1,j)) ...
                    + (eN(j)*p_k(i,j+1)+p_kp1(i,j-1))) ...
                    - g(i,j)*dX^2);

                % Calculating residual
                res(i,j) = (eE(i)*(p_k(i+1,j)-p_k(i,j)) ...
                    + eW(i)*(p_k(i-1,j)-p_k(i,j)) ...
                    + eN(j)*(p_k(i,j+1)-p_k(i,j)) ...
                    + (p_k(i,j-1)-p_k(i,j)))/dX^2 - g(i,j);
            end
        end
        % Calculating L inf norm of the residual
        B = reshape(res, [(I+1)*(J+1), 1]);
        L_inf_residual = norm(B, 'inf');

        % Setting the previous iteration equal to the new one
        p_k = p_kp1;
        iter = iter + 1;
        pois_total = pois_total + 1;
    end
    pois_iter(t) = iter; 
    iter = 0;

    %% Solving for u and v
    % Solving for u
    for j = 2:J
        for i = 2:I-1
            u(i,j) = u_star(i,j) - dT/(rho*dX)*(p_kp1(i+1,j)-p_kp1(i,j));
        end
    end
    % Solving for v
    for j = 2:J-1
        for i = 2:I
            v(i,j) = v_star(i,j) - dT/(rho*dY)*(p_kp1(i,j+1)-p_kp1(i,j));
        end
    end
end
%% Plots-------------------------------------------------------------------
fontSize = 12;
% Plot 1 - u velocity
figure("units","normalized","position",[0,0.33,0.3,0.3])
surf(xu(2:I-1),yu(2:J),u(2:I-1,2:J)')
xlabel('x')
ylabel('y')
set(gca,'fontsize',26)
title('u')
subtitle(['U = ', num2str(uLid), ', \gamma = ', num2str(gamma)], 'FontSize', fontSize);

% Plot 2 - v velocity
figure('units','normalized','position',[0,0.01,0.3,0.3])
surf(xv(2:I),yv(2:J-1),v(2:I,2:J-1)')
xlabel('x')
ylabel('y')
set(gca,'fontsize',26)
title('v')
subtitle(['U = ', num2str(uLid), ', \gamma = ', num2str(gamma)],'FontSize',  fontSize);

% Plot 3 - pressure
figure('units','normalized','position',[0.33,0.01,0.32,0.32])
surf(x_p(2:I),y_p(2:J),p_k(2:I,2:J)')
set(gca,'fontsize',26)
xlabel('x')
ylabel('y')
title('pressure')
subtitle(['U = ', num2str(uLid), ', \gamma = ', num2str(gamma)],'FontSize', fontSize);

% Plot 4 - u validation 
y_p_d=[1.0000 0.9766  0.9688 0.9609  0.9531 0.8516  0.7344 0.6172 0.5000 0.4531 0.2813 0.1719  0.1016 ...
0.0703 0.0625 0.0547 0.0000];%y coordinate
u_re100=[1.0000 0.8412 0.7887 0.7372 0.68717 0.2315 0.0033  -0.1364  -0.2058  -0.2109  -0.1566 ...
-0.1015  -0.0643  -0.04775  -0.0419  -0.0371 0.0000];% Re=100
figure
u_results = u(0.5/dX,:);
plot(y_p_d, u_re100)
hold on;
plot(y_p,u_results,"*")
hold off;
xlim([-0.2,1.2])
ylim([-0.4,1.2])
legend("Re=100 data", "Re=100 computation",'Location','northwest')
title('u vs Y comparison (at X=0.5)')
subtitle(['U = ', num2str(uLid), ', \gamma = ', num2str(gamma)],'FontSize', fontSize);
xlabel('y')
ylabel('u')

% Plot 5 - v validation 
% Experimental data
x_p_d=[1.0000 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5000 0.2344 0.2266 0.1563 0.0938 ...
0.0781 0.0703 0.0625 0.0000];
v_re100=[0.0000 -0.05906  -0.0739 -0.0886 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 ...
0.17507 0.16077 0.12317 0.1089 0.1009 0.0923 0.0000];
% My results
v_results = v(:,0.5/dY);
figure
hold on;
plot(x_p_d,v_re100)
plot(x_p,v_results,"*")
hold off;
xlim([-0.2,1.2])
ylim([-0.4,0.4])
legend("Re=100 data", "Re=100 computation",'Location','northwest')
title('v vs X comparison (at Y=0.5)')
subtitle(['U = ', num2str(uLid), ', \gamma = ', num2str(gamma)],'FontSize', fontSize);
xlabel('x')
ylabel('v')

% Plot 6 
figure;
plot(1:length(0:dT:tFinal), pois_iter)
title('Poisson Iterations vs. Time Iteration')
subtitle(['U = ', num2str(uLid), ', \gamma = ', num2str(gamma)]);
ylabel('Poisson Iterations')
xlabel('Time Iteration')

%% Function to generate U nodes
function [u,v,xu,yu,xv,yv,p,x_p,y_p,I,J] = generateNodes(dX, dY, L)
% Nodes in x
xu = 0:dX:L;
xv = -dX/2:dX:L+dX/2;

% Nodes in y
yu = -dY/2:dY:L+dY/2;
yv = 0:dY:L;

% Applying initial conditions
u = zeros(length(xu),length(yu));
v = zeros(length(xv),length(yv));

% Calculating indicies 
I = length(xu);
J = length(yv);

% Creating p points and grid
x_p = -dX/2:dX:L+dX/2;
y_p = -dY/2:dY:L+dY/2;
p = zeros(I+1,J+1);
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