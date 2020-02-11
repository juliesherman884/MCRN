%%% Constants and values given in the paper
% No solar forcing
tao1 = 122.63;
r = 0.62;
s = 1.97;
xi = 1.66e-3;
deltaw = 0.07;
zeta = 0.23;
tao2 = 30.37; 
eta = 86.9;
deltal = 5e-4;
ur0 = 10;
lambda0 = 0.75;
deltalambdaa = 2.25;
deltalambdas = 0.1;
zt = 50; %km

syms t

phi0dot = 0;
lambdadot = 0;
phi0 = 0;
ur = 0; % ur = ur0 - lambda*zt/2
lamda = lamnda0 + deltalambdaa*sin(2*pi*t)+ deltalambdas*sin(pi*t/11)^2; %annual and solar cycles


f = @(t,a) [-a(1)/tao1 - r*a(2)+ s*a(2)*a(3) - xi*phi0 + deltaw*phi0dot;
    -a(2)/tao1 + r*a(1)- s*a(1)*a(3) - zeta*phi0*a(3);
    -(a(3)-ur)/tao2 - eta*phi0*a(2) - deltal*lambdadot;];
xt0 = [10,10,10];
[tspan,a] = ode45(f,[0 100],xt0);     % Runge-Kutta 4th/5th order ODE solver
