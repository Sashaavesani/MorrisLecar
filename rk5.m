function [tout,yout] = rk5(odefun, tspan, y0, options)

% usando il tableau dell'esercizio implemento il metodo di Runge Kutta
% esplicito 
% odefun è la funzione
% tspan contiene solo due valori [t0, t0+tstar]
% y0 dato inziale 
% options contine options.InitalStep

c = [0;1/2;1/2;1;3/4];
b = [-1/2,7/3,7/3,13/6,-16/3];
a = [0,0,0,0,0;...
    1/2,0,0,0,0;...
    0,1/2,0,0,0;...
    0,0,1,0,0;...
    5/32,7/32,13/32,-1/32,0];

InitialStep = 0.01; % di defult, altrimenti sovrascritto

if (nargin == 4) % è stata passato l'input options
  if (isfield(options,'InitialStep'))
    InitialStep = options.InitialStep;
  end
end

k = InitialStep;
n = 1;                  % passo 1
tout(n) = tspan(1);
yout(:,n) = y0;
f = zeros(length(y0),5); % inizializzo le funzioni

while (tspan(2) - tout(n) > eps)
    k = min(k, tspan(2)-tout(n));
    f(:,1) = odefun(tout(n) + c(1)*k , yout(:,n));
    f(:,2) = odefun(tout(n) + c(2)*k, yout(:,n)+ a(2,1)*k*f(:,1));
    f(:,3) = odefun(tout(n) + c(3)*k, yout(:,n)+ a(3,2)*k*f(:,2));
    f(:,4) = odefun(tout(n) + c(4)*k, yout(:,n)+ a(4,3)*k*f(:,3));
    f(:,5) = odefun(tout(n) + c(5)*k, yout(:,n)+ k* (a(5,1)*f(:,1) + ...
        a(5,2)*f(:,2) + a(5,3)*f(:,3) + a(5,4)*f(:,4)));
    
    yout(:,n+1) = yout(:,n) + k*(b(1)*f(:,1)+b(2)*f(:,2)+b(3)*f(:,3)+b(4)*f(:,4)+b(5)*f(:,5));
    tout(n+1) = tout(n) + k;
    n = n+1;
end

tout = tout.';
yout = yout.';