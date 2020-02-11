%%%%%%%%%%%%%% Constants given in the paper %%%%%%%%%%%%%%
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
deltalambdas = 0.1; % No solar forcing if = 0
zt = 50; %km
g = 9.807; %m s-2
f0 = 1.26e-4; %per sec


%%%%%%%%%%%%%% Control Functions %%%%%%%%%%%%%%

lambda = @(t) lambda0 + deltalambdaa*sin(2*pi*t)+ deltalambdas*sin(pi*t/11)^2; %annual and solar cycles
lambdadot = @(t) 2*pi*deltalambdaa*cos(2*pi*t)+2*pi/11*deltalambdas*sin(pi*t/11)*cos(pi*t/11);
h = @(t) 0; %%%%%%%%%%%%%%%WHAT IS H
hdot = @(t) 0;

phi0 = @(t) g/f0*h(t);              
phi0dot = @(t) g/f0*hdot(t);
ur = @(t) ur0 - lambda(t)*zt/2;       % ur = ur0 - lambda*zt/2

%%%%%%%%%%%%%% Solve DE %%%%%%%%%%%%%%
f = @(t,a) [-a(1)/tao1 - r*a(2)+ s*a(2)*a(3) - xi*phi0(t) + deltaw*phi0dot(t);
    -a(2)/tao1 + r*a(1)- s*a(1)*a(3) - zeta*phi0(t)*a(3);
    -(a(3)-ur(t))/tao2 - eta*phi0(t)*a(2) - deltal*lambdadot(t);];
xt0 = [10,10,10];
[tspan,a] = ode45(f,[0 100],xt0);  
plot(tspan,a(:,3))
