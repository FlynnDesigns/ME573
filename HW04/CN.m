% Nathan Flynn
% Function for CN
function [T_CN,x] = CN(xStart, xEnd, kappa, deltaT, deltaX, time)
% Creating spacial nodes
t = 0:deltaT:time;
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