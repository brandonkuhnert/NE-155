function [ x,iter ] = Jacobi( A,b )
%Input: a matrix, the b vector, and number of iterations k, 
%value for x^(0)
%output: the x vector

% diagonal part of A and rest
D = diag(diag(A));
R = A - D;

% iteration matrix and offset
T = - inv(D) * R;
C = inv(D) * b;


% spectral radius condition
rho = max(abs(eig(T)));
if rho >= 1
    error('no convergence')
end

%number of iterations
iter=-6/(log10(rho));

% initial guess
x = randn(size(b));

% iteration
while norm(A * x - b) > 10e-6
    x = T * x + C;
end

