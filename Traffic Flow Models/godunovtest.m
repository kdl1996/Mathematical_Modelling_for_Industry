% Godunov scheme to solve burgers equation arising in the traffic flow
% model

clear all;

%umax=1 rhomax=1

xend = 15;
% x-axis size. 
tend = 5;
% t-axis size.
N = 15; 
dx = xend/N;
% Grid spacing 
dt = 1;

x = 0:dx:xend;
nt = floor(tend/dt); 
dt = tend / nt;
% initial traffic density values.
u0=[0.55,0.55,0.55,0.55,0.55,0.55,-1,0,0,0,0,0,0,0,0,0]; %for red light
%u0=[0,0,0,0,0,0,0,0,0,0,5.5,5.5,5.5,5.5,5.5,5.5];
%u0=[0.7866,0.93255,0.98469,0.99728,0.99965,0.99997,0,0,0,0,0,0,0,0,0,0]; %for green light

u = u0; 
unew = 0*u;
%Implementation of the Godunov Scheme.

    for i = 1 : nt,
        unew(2:end-1) =u(2:end-1)- dt/dx*(nf(u(2:end-1),u(3:end)) - nf(u(1:end-2),u(2:end-1))); % Different cases for the Godunov scheme is implemented in the function 'nf'
        unew(1) = u(1); 
        unew(end) = u(end);
       
        u = unew; 
       
        U(i,:) = u(:);
        fprintf('For t=%d \n',i);
        disp(U(i,:));
        
     
      
    end


U=[u0;U];
T=0:dt:tend;
%Plot of the solutions.
figure(1)
surf(x,T,U)
shading interp 
xlabel('x'),ylabel('t'), zlabel ('u(x,t)');
grid on 
colormap('parula');



                
