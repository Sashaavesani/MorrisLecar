% biforcation diagramm
clear all
clc

% problem setting and phase portrait
gca = 1;
gk = 2;
gl = 0.5;
Eca = 1;
Ek = -0.7;
El = -0.5;
phi = 1/3;
V1 = -0.01;
V2 = 0.15;
V3 = 0.1;
V4 = 0.145;

minf = @(V) 1/2 * (1 + tanh((V-V1)/(V2)));
winf = @(V) 1/2 * (1 + tanh((V-V3)/(V4)));
tau = @(V) (cosh((V-V3)/(2*V4)))^(-1);
Ica = @(V) gca*minf(V).*(V-Eca);
Il = @(V) gl*(V-El);


Iapp = linspace(0,0.4,100);

a = zeros(2,200);
b = zeros(1,200);

z = zeros(1,200);

counter = 0;
for i = Iapp
    f = @(Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+i;...
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
    e1 = e1+delta;
    
    if det(Jf(e1)) < Inf
        eig1 = eig(Jf(e1));
        if real(eig1(1)) > 0 && real(eig1(2)) > 0
            scatter(i,e1(1),'.r')
            counter = counter+1;
            b(counter) = i;
            a(:,counter) = e1;
            hold on
        end 
        if (real(eig1(1)) > 0 && real(eig1(2)) < 0) || (real(eig1(1)) < 0 && real(eig1(2)) > 0)
            scatter(i,e1(1),'.b')
            hold on
        end 
        if real(eig1(1)) < 0 && real(eig1(2)) < 0
            scatter(i,e1(1),'.k')
            hold on
        end
        if real(eig1(1)) == 0 && real(eig1(2)) == 0
            counter = counter + 1;
            z(counter) = i;
        end 
    end 

    y0 = [-0.16;0.02];
    e2 = y0;
    delta = -Jf(e2)\f(e2);
    while norm(delta,inf)>tol
        e2 = e2+delta;
        delta = -Jf(e2)\f(e2);
    end 
    e2 = e2+delta;
    if det(Jf(e2)) < Inf
        eig2 = eig(Jf(e2));
        if real(eig2(1)) > 0 && real(eig2(2)) > 0
            scatter(i,e2(1),'.r')
            counter = counter+1;
            b(counter) = i;
            a(:,counter) = e2;
            hold on
        end 
        if (real(eig2(1)) > 0 && real(eig2(2)) < 0) || (real(eig2(1)) < 0 && real(eig2(2)) > 0)
            scatter(i,e2(1),'.b')
            hold on
        end 
        if real(eig2(1)) < 0 && real(eig2(2)) < 0
            scatter(i,e2(1),'.k')
            hold on
        end 
        if real(eig2(1)) == 0 && real(eig2(2)) == 0
            counter = counter + 1;
            z(counter) = i;
        end 
    end 

    y0 = [0.03;0.27];
    e3 = y0;
    delta = -Jf(e3)\f(e3);
    while norm(delta,inf)>tol
        e3 = e3+delta;
        delta = -Jf(e3)\f(e3);
    end 
    e3 = e3+delta;
    if det(Jf(e3)) < Inf
        eig3 = eig(Jf(e3));
        if real(eig3(1)) > 0 && real(eig3(2)) > 0
            scatter(i,e3(1),'.r')
            counter = counter+1;
            b(counter) = i;
            a(:,counter) = e3;
            hold on
        end 
        if (real(eig3(1)) > 0 && real(eig3(2)) < 0) || (real(eig3(1)) < 0 && real(eig3(2)) > 0)
            scatter(i,e3(1),'.b')
            hold on
        end 
        if real(eig3(1)) < 0 && real(eig3(2)) < 0
            scatter(i,e3(1),'.k')
            hold on
        end 
        if real(eig3(1)) == 0 && real(eig3(2)) == 0
            counter = counter + 1;
            z(counter) = i;
        end 
    end 
    
end 

xlabel('I_{app}')
ylabel('V')
text(0.32,-0.1,'Instable Equilibia','Color','red')
text(0.32,-0.15,'Saddle Equilibria','Color','blue')
text(0.32,-0.20,'Stable Equilibria','Color','black')