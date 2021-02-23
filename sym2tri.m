function T = sym2tri(A)  % Function for Symmetric to Tridiagonal Conversion
n = size(A,1);

for k = 1:n-2             % Householder reduction to Hessenberg
    x = A(k+1:n, k);      
    sgn = sign(x(1,1));   
    e1 = zeros(n-k,1);    
    e1(1,1) = 1;          
    v = sgn*norm(x)*e1 + x;
    if norm(v) ~= 0
        v = v/norm(v);
    end
    A(k+1:n, k:n) = A(k+1:n, k:n) - 2*v*(v'*A(k+1:n, k:n));
    A(1:n, k+1:n) = A(1:n, k+1:n) - 2*(A(1:n, k+1:n)*v)*v';
end
T = A;     % Final Hessenberg form obtained is Tridiagonal and symmetric since A is symmetric
end