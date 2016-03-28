[A,b]=TridiagBuild(100);
P=inv(A);
x=P*b;
X=A\b;
plot(b,x,'ro',b,X,'b')
legend('exact solution','estimated solution')
ylabel('Solutions (x vectors)')
xlabel('b vector')