n=10
x=rand(n-1,1)
x=[0;sort(x);1]
A=zeros(n-1)
A=A+diag(1./(x(3:end)-x(2:end-1))+1./(x(2:end-1)-x(1:end-2)))
A=A+diag(1./(x(3:end-1)-x(2:end-2)),1)
A=A+diag(1./(x(3:end-1)-x(2:end-2)),-1)
A=A-2*diag(1./(x(3:end-1)-x(2:end-2)),-1)
A=A-2*diag(1./(x(3:end-1)-x(2:end-2)),1)
b=2*0.5*(x(3:end)-x(1:end-2))
u=A\b

% Numerical solution
y(:,1)=[0;u;0]
% Actual soluion
y(:,2)=-x.*x+x  

plot(x,y(:,1),'-o')
hold on
plot(x,y(:,2),'-x')
title('Comparison of actual and numerical solution')