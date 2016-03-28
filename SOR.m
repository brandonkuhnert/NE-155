function [ x ] = SOR( A,b,k,x,w )
%Input: a matrix, the b vector, and number of iterations k, 
%value for vector x^(0), w control value
%output: the x vector

if x==0;
    x=zeros(length(A),1);
    D=diag(diag(A));
    [L,U]=lu(A);
    F=inv(D+w.*L);
    
    for i=1:k;
        x=F*((1-w)*D-w.*U)*x+w.*F*b;
    end
end
end