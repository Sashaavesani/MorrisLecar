% andamento orbite con PP006


clear all
clc

%% problem setting and phase portrait
gca = 1;
gk = 2;
gl = 0.5;
Eca = 1;
Ek = -0.7;
El = -0.5;
phi = 1/3;
Iapp = 0.06;
V1 = -0.01;
V2 = 0.15;
V3 = 0.1;
V4 = 0.145;

minf = @(V) 1/2 * (1 + tanh((V-V1)/(V2)));
winf = @(V) 1/2 * (1 + tanh((V-V3)/(V4)));
tau = @(V) (cosh((V-V3)/(2*V4)))^(-1);
Ica = @(V) gca*minf(V).*(V-Eca);
Il = @(V) gl*(V-El);

f = @(t,Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
    phi*(winf(Y(1))-Y(2))/tau(Y(1))];

y1 = linspace(-0.6,0.3,20);
y2 = linspace(-0.05,0.40,20);
[x,y] = meshgrid(y1,y2);

u = zeros(size(x));
v = zeros(size(x));

t=0; 
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
u = u./(sqrt(u.^2+v.^2));
v = v./(sqrt(u.^2+v.^2));
q = quiver(x,y,u,v,0.5,'k', 'LineWidth',0.4);
figure(gcf)
xlabel('V')
ylabel('W')

hold on
W_nullocline = winf(y1);
plot(y1,W_nullocline,'c--', 'LineWidth',1);

hold on
V_nullocline = (Iapp-Ica(y1)-Il(y1))./(gk.*(y1-Ek));
plot(y1,V_nullocline,'r--','color', [1 .5 0], 'LineWidth',1);

xlabel('V')
ylabel('W')
xlim([-0.6,0.3])
ylim([-0.05,0.40])

%% fixed point with Newton
f = @(Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
    phi*(winf(Y(1))-Y(2))/tau(Y(1))];

h1 = @(x) tanh((x-V1)./(V2));
h2 = @(x) (1./V2) * (1-(tanh((x-V1)./(V2))).^2);
h3 = @(x) tanh((x-V3)./(V4));
h4 = @(x) (1./V4) * (1-(tanh((x-V3)./(V4))).^2);

c = @(x) cosh( (x-V3)./(2*V4) );
s = @(x) (1./(2*V4)) * sinh( (x-V3)./(2*V4) );

J1 = @(Y) -gca*(1/2*h2(Y(1)))*(Y(1)-Eca)-gca*(1/2+1/2*h1(Y(1)))-gk*Y(2)-gl;
J2 = @(Y) -gk*(Y(1)-Ek);
J3 = @(Y) phi.*(1/2*h4(Y(1))*c(Y(1))+(1/2+1/2*h3(Y(1))-Y(2))*s(Y(1)));
J4 = @(Y) -phi./tau(Y(1));
Jf = @(Y) [J1(Y),J2(Y);J3(Y),J4(Y)];


tol = 1e-10;
y0 = [-0.3535;0.002];
e1 = y0;
delta = -Jf(e1)\f(e1);
while norm(delta,inf)>tol
    e1 = e1+delta;
    delta = -Jf(e1)\f(e1);
end 
e1 = e1+delta
hold on
scatter(e1(1),e1(2),'g','filled', 'LineWidth',1.5)
eigvalues_e1 = eig(Jf(e1))

y0 = [-0.16;0.02];
e2 = y0;
delta = -Jf(e2)\f(e2);
while norm(delta,inf)>tol
    e2 = e2+delta;
    delta = -Jf(e2)\f(e2);
end 
e2 = e2+delta
hold on
scatter(e2(1),e2(2),'y','filled', 'LineWidth',1.5)
eigvalues_e2 = eig(Jf(e2))

y0 = [0.03;0.27];
e3 = y0;
delta = -Jf(e3)\f(e3);
while norm(delta,inf)>tol
    e3 = e3+delta;
    delta = -Jf(e3)\f(e3);
end 
e3 = e3+delta

hold on
scatter(e3(1),e3(2),'m','filled', 'LineWidth',1.5)

eigvalues_e3 = eig(Jf(e3))

%% Draw some orbits
f = @(t,Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
    phi*(winf(Y(1))-Y(2))/tau(Y(1))];

y0 = [-0.45;0.32];
[tout1,yout1] = ode45(f,[0,50],y0);
hold on
scatter(-0.45,0.32,'r', 'filled')
hold on
plot(yout1(:,1),yout1(:,2), 'r',1.2);
% for j = 1:length(yout1)-1
%     plot(yout1(j:j+1,1),yout1(j:j+1,2),'color', [1 .5 0], 'LineWidth',2.5)
%     hold on
%     drawnow
%     pause(.05)
% end

y0 = [-0.12;0.12];
[tout2,yout2] = ode45(f,[0,50],y0);
hold on
scatter(-0.12,0.12,'b', 'filled')
hold on
plot(yout2(:,1),yout2(:,2),'b',1.2)
% for j = 1:length(yout2)-1
%     plot(yout2(j:j+1,1),yout2(j:j+1,2),'m', 'LineWidth',2.5)
%     hold on
%     drawnow
%     pause(.05)
% end

y0 = [0.005;0.23];
[tout3,yout3] = ode45(f,[0,50],y0);
hold on
scatter(0.005,0.23,'k', 'filled')
hold on
plot(yout3(:,1),yout3(:,2), 'k',1.2)
% for j = 1:length(yout3)-1
%     plot(yout3(j:j+1,1),yout3(j:j+1,2),'k', 'LineWidth',2.5)
%     hold on
%     drawnow
%     pause(.05)
% end

legend('vectorial field','W nullocline', 'V nullocline', 'equilibrium 1', 'equilibrium 2', 'equilibrium 3', 'Location', 'southeast');

%% tV 
figure
plot(tout1, yout1(:,1),'r',  'LineWidth',1.5)
hold on 
plot(tout2, yout2(:,1),'b', 'LineWidth',1.5)
hold on 
plot(tout3, yout3(:,1),'k', 'LineWidth',1.5);
title('Time-Potential plot of three different orbits')
legend('V1', 'V2', 'V3')
xlabel('t')
ylabel('V (potential)')

%% tW 
figure
plot(tout1, yout1(:,2),'r',  'LineWidth',1.5)
hold on 
plot(tout2, yout2(:,2),'b', 'LineWidth',1.5)
hold on 
plot(tout3, yout3(:,2),'k', 'LineWidth',1.5);
title('Time-W plot of three different orbits')
legend('W1', 'W2', 'W3')
xlabel('t')
ylabel('W')