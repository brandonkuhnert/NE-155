function [ x ] = GaussSeidel( A,b,k,x )
%Input: a matrix, the b vector, and number of iterations k, 
%value for vector x^(0)
%output: the x vector

if x==0;
    x=zeros(length(A),1);
    D=diag(diag(A));
    [L,U]=lu(A);
    F=inv(D+L);
    
    for i=1:k;
        x=-U*F*x+F*b;
    end
end


end

