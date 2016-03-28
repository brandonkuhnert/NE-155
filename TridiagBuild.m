function [ A,b ] = TridiagBuild( n )
%Brandon Kuhnert

j=[1:n-2];
b=[0,j,n-1]';
m=ones(n,1)+ones(n,1);
d=-1.*ones(n-1,1);
A=diag(m)+diag(d,-1)+diag(d,1);
end

