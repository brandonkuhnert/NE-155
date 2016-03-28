function [ x ] = Jacobi( A,b,k,x )
%Input: a matrix, the b vector, and number of iterations k, 
%value for x^(0)
%output: the x vector

if x==0;
    x=zeros(length(A),1);
    D=diag(diag(A));
    F=inv(D);
    
    for i=1:k;
        x=(D*F-F*A)*x+F*b;
    end
end

end

