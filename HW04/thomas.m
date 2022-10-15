function x = thomas(A,b,N)
x = zeros(N,1);
% Transforming into an upper triangular form
for k = 1:(N-1)
    i = k+1;
    l_ik = A(i,k)/A(k,k);
    for j = k:(k+1)
        A(i,j)= A(i,j)-l_ik*A(k,j);
    end
    b(i) = b(i)-l_ik*b(k);
end
% Applying backwards substitution 
for k = N:-1:1
    x(k) = b(k);
    for j = k+1:min(N,k+1)
        x(k) = x(k) - A(k,j)*x(j);
    end
    x(k) = x(k)/A(k,k);
end 
