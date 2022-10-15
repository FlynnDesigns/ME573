% Nathan Flynn
% 10/02/2022
% ME573 
% HW04

% Problem01
clc; close all; clear;

% Constants
xStart = -3;
xEnd = 3;
kappa = 5 * 10^-3;
deltaX = 0.1;
deltaT = 0.1;
time = 100;
alpha =  (kappa*deltaT)/(deltaX^2);

% FTCS
[T_FTCS, x_FTCS] = FTCS(xStart, xEnd, kappa, deltaT, deltaX, time);

% CN
[T_CN, x_CN] = CN(xStart, xEnd, kappa, deltaT, deltaX, time);

% BTCS
[T_BTCS, x_BTCS] = BTCS(xStart, xEnd, kappa, deltaT, deltaX, time);

% Analytical solution
x = xStart : deltaX : xEnd;
t = 0:(time/deltaT);
Texact = zeros(1, length(x));

for i = 1:(length(x)) % Sweeps through spacial nodes
    Texact(i) = (erf((1-x(i))/(2*sqrt(kappa*(time)))) - erf(-(x(i)+1)/(2*sqrt(kappa*time)))); 
end

Linf_FTCS = norm(T_FTCS - Texact, Inf);
Linf_BTCS = norm(T_BTCS - Texact, Inf);
Linf_CN = norm(T_CN - Texact, Inf);

fprintf('deltaX = %.2d, deltaT = %.2d, and time = %.0d \n',deltaX, deltaT, time);
fprintf('Linf_{FTCS} = %d \n', Linf_FTCS);
fprintf('Linf_{BTCS} = %d \n', Linf_BTCS);
fprintf('Linf_{CN} = %d \n', Linf_CN);


% Part b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting
figure('units','normalized','position',[0.55 0.1 0.45 0.45]);
plot(x,Texact,'-o',x,T_FTCS,'-d',x,T_BTCS,'--',x,T_CN,'-+');
ylim([0,3]);
ax = gca;
set(gca,'fontsize',26);
title(sprintf("\\Delta x = %.2d , \\Delta t = %.2d , time = %.0d" ,deltaX, deltaT, time));
ax.TitleFontSizeMultiplier = 0.5;
legend('Exact','FTCS','BTCS','Crank-Nicolson');
xlabel('x')

% Checking stability for FTCS
if alpha <= 1/2
    disp("FTCS is stable");
else
    disp("FTCS is unstable");
end

% BTCS is unconditinally stable 

% CN is unconditionally stable

%% Function for FTCS
function [T_FTCS,x] = FTCS(xStart, xEnd, kappa, deltaT, deltaX, time)
% Creating spacial nodes
x = xStart:deltaX:xEnd;
alpha =  (kappa*deltaT)/(deltaX^2);

% Creating vectors
fn = zeros(1,length(x));
fnp1 = zeros(1,length(x));

% Applying initial conditions
for i = 1:length(x)
    if abs(x(i)) <= 1
        fn(i) = 2;
    else
        fn(i) = 0;
    end
end

% Solving using FTCS
for j = 1:(length(0:deltaT:time)-1) % Sweeps through time nodes
    for i = 2:(length(x)-1) % Sweeps through spacial nodes
        fnp1(i) = fn(i) + alpha*(fnp1(i+1) - 2*fn(i) +fn(i-1));
    end
    fn = fnp1;
end
T_FTCS = fn;
end

%% Function for CN
function [T_CN,x] = CN(xStart, xEnd, kappa, deltaT, deltaX, time)
% Creating spacial nodes
x = xStart:deltaX:xEnd;
N = length(x);
alpha =  (kappa*deltaT)/(deltaX^2);
e = ones(N-2,1);
A = spdiags([-alpha*e 2*(alpha+1)*e -alpha*e],-1:1,N-2,N-2);

% Creating vectors
T_CN = zeros(1,length(x));

% Applying initial conditions
for i = 1:N
    if abs(x(i)) <= 1
        T_CN((i)) = 2;
    else
        T_CN((i)) = 0;
    end
end

for t = 0:deltaT:time
    b = alpha * T_CN(1:N-2) + 2*(1-alpha) * T_CN(2:N-1) + alpha * T_CN(3:N);
    b(1) = b(1) + alpha * T_CN(1);
    b(N-2) = b(N-2) + alpha * T_CN(N);
    z = thomas(A,b,N-2);
    T_CN(2:N-1) = z(1:N-2);
end
    
end

%% Function for BTCS
function [T_BTCS,x] = BTCS(xStart, xEnd, kappa, deltaT, deltaX, time)
% Creating spacial nodes
x = xStart:deltaX:xEnd;
N = length(x);
alpha =  (kappa*deltaT)/(deltaX^2);
e = ones(N,1);
A = spdiags([-alpha*e (2*alpha+1)*e -alpha*e],-1:1,N-2,N-2);

% Creating vectors
T_BTCS = zeros(1,length(x));
b = zeros(1,length(x));

% Applying initial conditions
for i = 1:N
    if abs(x(i)) <= 1
        T_BTCS((i)) = 2;
    else
        T_BTCS((i)) = 0;
    end
end

for t = 0:deltaT:time
    for i = 2:N-2
        b(i) = T_BTCS(i+1);
    end
    b(1) = T_BTCS(2) + alpha * T_BTCS(1);
    b(N) = T_BTCS(N-1) + alpha * T_BTCS(N);
    z = thomas(A,b,N-2);
    T_BTCS(2:N-1) = z(1:N-2);
end
    
end