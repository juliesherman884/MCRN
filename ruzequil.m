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



%%%%%%%%%%%%%% Equilibrium %%%%%%%%%%%%%%
hdot = 0;
lambdadot = 0;
%Fixed Lambda
lambdae = 1;
lambda = lambdae;
h = [0:10:250];
x = zeros(length(h),3);

for i = 1:length(h)
    ht = h(i);
    phi0 = g/f0*ht;
    fun = @(a) [-a(1)/tao1 - r*a(2)+ s*a(2)*a(3) - xi*phi0;
        -a(2)/tao1 + r*a(1)- s*a(1)*a(3) - zeta*phi0*a(3);
        -(a(3)-ur)/tao2 - eta*phi0*a(2);];
    x0 = [0,0,30];
    x(i,:) = fsolve(fun,x0);
end
plot(h,x(:,3))