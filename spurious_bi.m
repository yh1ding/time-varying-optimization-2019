%% TAC paper code
% Code for numerical example in the paper 'On the Absence of Spurious Local 
% Trajectories in Time-varying Nonconvex Optimization'
% Author: Salar Fattahi, Cedric Josz, Yuhao Ding, Reza Mohammadi, Javad Lavaei, and Somayeh Sojoudi

P = 1;
T = P*2*pi;
a = .4; 
b = 10;
[X,Y] = meshgrid(-4:.1:4,0:.1*P:T);
x = X - b*sin(Y)+(b-1)*sin(Y);
Z = (1/4).*x.^4 + 1/8.*x.^3 + -2*x.^2  -3/2*x + 8;
Z = min(Z,9);

steps = 10000*P;
dt = T/steps;
x = zeros(steps+1,1);
y = zeros(steps+1,1);
z = zeros(steps+1,1);
x(1) = -2;
% x(1) = 2;
for k = 1:steps
    t = (k-1)*dt;
    x(k+1) = x(k) - dt/a * ( (x(k)-b*sin(t))^3 + 3/8*(x(k)-b*sin(t))^2 - 4*(x(k)-b*sin(t)) -3/2 );
    y(k) = t;
    z(k) = 1/4*(x(k)-b*sin(t))^4 + 1/8*(x(k)-b*sin(t))^3 - 2*(x(k)-b*sin(t))^2 -3/2*(x(k)-b*sin(t)) + 8;
end
k = steps+1; t = (k-1)*dt;
y(k) = t;
z(k) = 1/4*(x(k)-b*sin(t))^4 + 1/8*(x(k)-b*sin(t))^3 - 2*(x(k)-b*sin(t))^2 -3/2*(x(k)-b*sin(t)) + 8;
z = min(z,9);

for k=1:steps
    t = (k-1)*dt;    
    x(k) = x(k) - (b-1)*sin(t);
end

figure
surf(X,Y,Z)
title('Non-spurious Trajectory','FontWeight','bold')
xlabel('x: variable','FontWeight','bold')
ylabel('t: time','FontWeight','bold') 
zlabel('f(x,t)','FontWeight','bold') 
set(gca,'FontSize',15);
hold on
plot3(x,y,z,'.','MarkerSize',15,'Color','red')

% figure
% for k=1:steps
%     t = (k-1)*dt;    
%     x(k) = x(k) + (b-1)*sin(t);
% end
% for k = 1:size(z,1)
%     t = (k-1)*dt;
%     z(k,1) = z(k)+ (x(k)-(2+b*sin(t)))^2;% + 1/10*(x(k)-(2))^2;
% end
% plot(1:size(z,1),z,'Linewidth',3)
% xlabel('t: time','FontWeight','bold')
% ylabel('f(x,t)+||x-2-b sin(t)||^2','FontWeight','bold') 
