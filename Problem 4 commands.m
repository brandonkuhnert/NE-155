n=5
b=ones(n,1)
m=ones(n,1).*4
d=-1.*ones(n-1,1)
A=diag(m)+diag(d,-1)+diag(d,1)
