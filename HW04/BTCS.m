% Nathan Flynn
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