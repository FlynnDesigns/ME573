% Nathan Flynn
% Function for FTCS
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