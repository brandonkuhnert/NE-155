function [ x,iter,rho ] = SORrel( n,w,error )
%Input: a matrix, the b vector, and number of iterations k, 
%value for vector x^(0), w control value
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
T = inv(D+w.*L);
E = ((1-w).*D-w.*U);
C = T * (w.*b);


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
while norm(x-oldx)/norm(x) > error
    iter=iter+1;
    oldx=x;
    x=T*E*oldx+C;
end
end