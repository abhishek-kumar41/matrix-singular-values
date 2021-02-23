clc 
clear all

A = [3, 3, 4;
    3, 7, 6;
    4, 6, 10];

n = size(A);

S = A'*A;  %Converting given matrix to symmetric matrix by multiplying transpose of A with A

T = sym2tri(S);  % Symmetric to Tridiagonal Conversion

T1 = zeros(n);   
T2 = T;          % Values of T used for iteration

while norm(T2 - T1) >1e-20     
    T1 = T2;
    [Q,R] = Givens_rotation(T1);  %Computing QR factorisation of Tridiagonal matrix at each iteration
    T2 = R*Q;                     % T2 is similar to T1 as T2 = Q'*T1*Q
end

eigen_values_AtA = diag(T2);     %Eigenvalues of transpose(A)*A at diagonal of T2 since T2 is diagonal matrix
singular_values_A = sqrt(eigen_values_AtA); % Singular Values is square root of eigenvalues of transpose(A)*A
disp('Singular Values of Matrix is:')
disp(singular_values_A)