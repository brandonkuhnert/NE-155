function [ I ] = CompSimp38( a,b,n )
%Brandon Kuhnert
fun= @(x)x/((x^2-4)^(1/2));
h=(b-a)/n;
K = feval(fun,a);
for i=1:3:n-1
    x(i)=a+h*i;
    K=K+3*feval(fun,x(i));
end
for i=2:3:n-1
    x(i)=a+h*i;
    K=K+3*feval(fun,x(i));
end
for i=3:3:n-1
    x(i)=a+h*i;
    K=K+2*feval(fun,x(i));
end
K=K+feval(fun,b);
I=3*h*K/8;

end

