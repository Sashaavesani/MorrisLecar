
clear all
clc

%% problem setting
gca = 1;
gk = 2;
gl = 0.5;
Eca = 1;
Ek = -0.7;
El = -0.5;
phi = 1/3;
Iapp = 0.1;
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

y1 = linspace(-0.5,0.3,20);
y2 = linspace(-0.05,0.5,20);
[x,y] = meshgrid(y1,y2);

u = zeros(size(x));
v = zeros(size(x));

%% plot vectorial field through "quiver"
t=0; 
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
u = u./(sqrt(u.^2+v.^2));
v = v./(sqrt(u.^2+v.^2));
q = quiver(x,y,u,v,0.5,'color',[0.6 0.6 0.6], 'LineWidth',0.4);
figure(gcf)
xlabel('V')
ylabel('W')

hold on
W_nullocline = winf(y1);
plot(y1,W_nullocline,'c--', 'LineWidth',1.3);

hold on
V_nullocline = (Iapp-Ica(y1)-Il(y1))./(gk.*(y1-Ek));
plot(y1,V_nullocline,'r--','color', [1 .5 0], 'LineWidth',1.3);

xlabel('V')
ylabel('W')
xlim([-0.5,0.3])
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
y0 = [0.05;0.03];
e1 = y0;
delta = -Jf(e1)\f(e1);
while norm(delta,inf)>tol
    e1 = e1+delta;
    delta = -Jf(e1)\f(e1);
end 
e1 = e1+delta
hold on
scatter(e1(1),e1(2),'g','filled', 'LineWidth',1.2)
eigvalues_e1 = eig(Jf(e1))

%% Draw some orbits
f = @(t,Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
    phi*(winf(Y(1))-Y(2))/tau(Y(1))];

y0 = [0.005;0.29];
options.InitialStep = 0.001;
[tout1,yout1] = rk5(f, [0,100], y0, options);
hold on
scatter(0.005,0.29,'r','filled')
hold on
plot(yout1(:,1),yout1(:,2),'r',1.8);
hold on
plot(yout1(5000:200:6000,1),yout1(5000:200:6000,2),'r',1.3)
xlim([-0.5,0.3])
ylim([-0.05,0.5])

y0 = [-0.12;0.4];
options.InitialStep = 0.001;
[tout2,yout2] = rk5(f, [0,100], y0, options);
hold on
scatter(-0.12,0.4,'k','filled')
hold on
plot(yout2(:,1),yout2(:,2),'k',1.2);
hold on
plot(yout2(1:100:1000,1),yout2(1:100:1000,2),'k',1);
xlim([-0.5,0.3])
ylim([-0.05,0.5])

legend('vectorial field', 'W nullocline', 'V nullocline', 'equilibrium', 'Location', 'northwest')

%% PLOT t-V
figure
plot(tout1, yout1(:,1),'r','LineWidth',1.2);
hold on 
plot(tout2, yout2(:,1),'k','LineWidth',1.2);
legend('V1','V2')
title('Time-Potential plot of the orbit')
xlabel('t')
ylabel('V (potential)')

%% PLOT t-W
figure
plot(tout1, yout1(:,2),'r','LineWidth',1.2);
hold on 
plot(tout2, yout2(:,2),'k','LineWidth',1.2);
legend('W1','W2')
title('Time-W plot of the orbit')
xlabel('t')
ylabel('W')
