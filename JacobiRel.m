function [ x,iter ] = JacobiRel( n,error )
%Input: a matrix, the b vector, and number of iterations k, 
%value for x^(0)
%output: the x vector

b=ones(n,1).*100;
m=ones(n,1).*4;
d=-1.*ones(n-1,1);
A=diag(m)+diag(d,-1)+diag(d,1);

% diagonal part of A and rest
D = diag(diag(A));
R = A - D;

% iteration matrix and offset
T = - inv(D) * R;
C = inv(D) * b;


% spectral radius condition
rho = max(abs(eig(T)));
if rho >= 1
    error('no convergence');
end

%number of iterations
iter=-6/(log10(rho));
oldx=zeros(n,1);
iter=0;

% initial guess
x = randn(size(b));

% iteration
while norm(x-oldx)/norm(x) > error;
    iter=iter +1;
    oldx=x;
    x=T*x+C;
end