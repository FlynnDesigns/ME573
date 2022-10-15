% Nathan Flynn
% 10/08/2022
% ME 573 - HW 5
clc; clear; close all;
%% Part a - create a surface plot of |f_{series}(x, y, 0) - f_{init}(x,y)| ** value should be around 10^-6         
% Initial variables
kappa = 0.1;
deltaX = 0.05;
deltaY = 0.05;
nSeries = 100;

% Computing other initial variables
deltaT = (deltaX^2 / kappa) / 4;
tf = 20 * deltaT;
tInit = 0;
x = 0:deltaX:1;
y = 0:deltaY:1;
t = 0:deltaT:tf;
nX = length(x);
nY = length(y);
[X, Y] = meshgrid(x, y);
X = X';
Y = Y';

% Initializing arrays and applying initial conditions
f_exact = zeros(nX, nY);
f_init = X.*(1-X.^5) .* Y.*(1-Y);

% Main loop for cacluating f_exact
for n = 1:nSeries 
    for m = 1:nSeries
            top = 120*(-(n^4)*(pi^4)*(-1)^n + 12*(n^2)*(pi^2)*(-1)^n + 24 + 24*(-1)^(1 + n)) * (-2 + 2 * (-1)^m);
            bottom = (n^7)*(pi^10)*(m^3);
            coeff = sin(n*pi*X).*sin(m*pi*Y)*exp(-(n^2 + m^2) * pi^2 * kappa * tInit);
            f_exact = f_exact - (top * coeff) / bottom;
    end
end

% Plotting figures
figure('units', 'normalized','position', [0 0.01 .4 .4]);
surf(X,Y, abs(f_exact - f_init));
set(gca, 'fontsize', 18)
title('|f_{series}(x,y,0)-f_{init}(x,y)|')

%% Part b
alpha = (kappa * deltaT) / (deltaX^2);
ftcs_0 = f_init;
ftcs_1 = zeros(nX, nY);
            
for j = 1:(tf/deltaT)
     ftcs_1(2:nX-1, 2:nY-1) = ftcs_0(2:nX-1, 2:nY-1) + alpha * (ftcs_0(2:nX-1, 3:nY) - 2*ftcs_0(2:nX-1, 2:nY-1) + ftcs_0(2:nX-1, 1:nY-2)) + alpha * (ftcs_0(3:nX, 2:nY-1) - 2*ftcs_0(2:nX-1, 2:nY-1) + ftcs_0(1:nX-2, 2:nY-1));
     ftcs_0 = ftcs_1;
end

% Initializing arrays and applying initial conditions
f_series = zeros(nX, nY);

% Main loop for cacluating f_exact at t final
for n = 1:nSeries 
    for m = 1:nSeries
        % Solving the exact equation
        top = 120*(-(n^4)*(pi^4)*(-1)^n + 12*(n^2)*(pi^2)*(-1)^n + 24 + 24*(-1)^(1 + n)) * (-2 + 2 * (-1)^m);
        bottom = (n^7)*(pi^10)*(m^3);
        coeff = sin(n*pi*X).*sin(m*pi*Y)*exp(-(n^2 + m^2) * pi^2 * kappa * tf);
        f_series = f_series - (top * coeff) / bottom;
    end
end

% Plotting figures
figure('units', 'normalized','position', [0 0.01 .4 .4]);
surf(X,Y, abs(ftcs_1));
set(gca, 'fontsize', 18)
title('f_{FTCS}(x,y,t_{final})') % FTCS

%% Part c
% Plotting figures
figure('units', 'normalized','position', [0 0.01 .4 .4]);
surf(X,Y, abs(ftcs_1 - f_series));
set(gca, 'fontsize', 18)
title('|f_{FTCS}(x,y,t_{final})-f_{series}(x,y,t_{final})|') % FTCS vs exact

%% Part d
% Linf_FTCS = norm(T_FTCS - Texact, Inf);
B = reshape(ftcs_1 - f_series, [nX * nY,1]);
fprintf("L_inf norm of error = %f\n", norm(B,"inf"));