function [Q,R] = Givens_rotation(A)
n = size(A,1);

Q = eye(n);

for i = 1:n-1
    if A(i+1,i)~=0     %When desired term is already zero, no need to go through iteration
        G = eye(n);
        a = A(i,i);
        b = A(i+1,i);         %element in A which had to be made zero
        c = sqrt(a^2 + b^2);
        cos_theta = a/c;
        sin_theta = b/c;
        G(i,i) = cos_theta;         %Forming G' matrix of size 2x2 and rest elements of G at diagonal are 1
        G(i+1,i+1) = cos_theta;
        G(i,i+1) = sin_theta;
        G(i+1,i) = -sin_theta;
        A = G*A;                    %This will give final R
        Q = Q*G';                   % This will give final Q
    end
end
R = A;
end