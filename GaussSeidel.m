function [ x,iter ] = GaussSeidel( n )
%Input: a matrix, the b vector, and number of iterations k, 
%value for vector x^(0)
%output: the x vector
b=ones(n,1).*100;
m=ones(n,1).*4;
d=-1.*ones(n-1,1);
A=diag(m)+diag(d,-1)+diag(d,1);
x=zeros(n,1);

% diagonal part of A and LU decomposition
D = diag(diag(A));
[L,U]=lu(A);
U=U-diag(diag(U));
L=A-D-U;

% iteration matrix and offset
T =  inv(D+L) * -U;
C = inv(D+L) * b;


% spectral radius condition
rho = max(abs(eig(T)));
if rho >= 1
    error('no convergence')
end

%number of iterations
iter=-6/(log10(rho));

% initial guess
oldx = zeros(size(b));
oldx=ones(n,1).*100;
iter=0;

% iteration
while norm(x-oldx) > 10e-6
    iter=iter+1;
    oldx=x;
    x=T*oldx+C;
end

