clc;                                               
clear all;
h=1;                                             % step size
t = 0:h:20;                                         
V = zeros(1,length(t))
C = zeros(1,length(t))
F = zeros(1,length(t))
m = zeros(1,length(t))

% Model parameters

beta=0.17
a=1
gamma=10
mu=0.5
K=10000
sigma=10
eta=0.12
tau=0.5
Cstar=10
alpha=1
p=0.8
qk=1
tshift=tau/h



% Initial conditions.
V(1)= 0.000001;
F(1)= 10;
C(1)= 10;
m(1)= 0.1;

%Vs=circshift(V,[0,-tshift])
%Fs=circshift(F,[0,-tshift])

% The functions on the RHS of ode
fV = @(t,v,c,f,m) alpha*v-p*f*v
fF = @(t,v,c,f,m) beta*c-gamma*p*f*v-a*f  
fC = @(t,v,c,f,m) -mu*(c-Cstar)+qk*v*f 
fm = @(t,v,c,f,m) sigma*v-eta*m

% 4th order RK Methos implementation
for i=1:(length(t)-1)                              
    k1v = fV(t(i),V(i),C(i),F(i),m(i));
    k1f = fF(t(i),V(i),C(i),F(i),m(i));
    k1c = fC(t(i),V(i),C(i),F(i),m(i));
    k1m = fm(t(i),V(i),C(i),F(i),m(i));
    
    
    k2v = fV(t(i)+0.5*h,V(i)+0.5*h*k1v,C(i)+0.5*h*k1c,F(i)+0.5*h*k1f,m(i)+0.5*h*k1m);
    k2f = fF(t(i)+0.5*h,V(i)+0.5*h*k1v,C(i)+0.5*h*k1c,F(i)+0.5*h*k1f,m(i)+0.5*h*k1m);
    k2c = fC(t(i)+0.5*h,V(i)+0.5*h*k1v,C(i)+0.5*h*k1c,F(i)+0.5*h*k1f,m(i)+0.5*h*k1m);
    k2m = fm(t(i)+0.5*h,V(i)+0.5*h*k1v,C(i)+0.5*h*k1c,F(i)+0.5*h*k1f,m(i)+0.5*h*k1m);
    
    k3v = fV(t(i)+0.5*h,V(i)+0.5*h*k2v,C(i)+0.5*h*k2c,F(i)+0.5*h*k2f,m(i)+0.5*h*k2m);
    k3f = fF(t(i)+0.5*h,V(i)+0.5*h*k2v,C(i)+0.5*h*k2c,F(i)+0.5*h*k2f,m(i)+0.5*h*k2m);
    k3c = fC(t(i)+0.5*h,V(i)+0.5*h*k2v,C(i)+0.5*h*k2c,F(i)+0.5*h*k2f,m(i)+0.5*h*k2m);
    k3m = fm(t(i)+0.5*h,V(i)+0.5*h*k2v,C(i)+0.5*h*k2c,F(i)+0.5*h*k2f,m(i)+0.5*h*k2m);
    
    k4v = fV(t(i)+h,V(i)+h*k3v,C(i)+h*k3c,F(i)+h*k3f,m(i)+h*k3m);
    k4f = fV(t(i)+h,V(i)+h*k3v,C(i)+h*k3c,F(i)+h*k3f,m(i)+h*k3m);
    k4c = fV(t(i)+h,V(i)+h*k3v,C(i)+h*k3c,F(i)+h*k3f,m(i)+h*k3m);
    k4m = fV(t(i)+h,V(i)+h*k3v,C(i)+h*k3c,F(i)+h*k3f,m(i)+h*k3m);
    
    %Updating step
    V(i+1) = V(i) + (1/6)*(k1v+2*k2v+2*k3v+k4v)*h;
    F(i+1) = F(i) + (1/6)*(k1f+2*k2f+2*k3f+k4f)*h;
    C(i+1) = C(i) + (1/6)*(k1c+2*k2c+2*k3c+k4c)*h;
    m(i+1) = m(i) + (1/6)*(k1m+2*k2m+2*k3m+k4m)*h;
end
plot(t,V)
title('Plot for concentration of antigens')