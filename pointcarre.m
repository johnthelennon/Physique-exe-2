data = load('poincarre.out');
tet=data(:,2);
theta = wrapToPi(tet);
%theta = mod(tet,2*pi);
theta_point = data(:,3);

% r = 0.14
% L = 0.18
% Omega =7.382411530116701
% t= data(:,1)
% x = r*cos(Omega*t) + L*sin(theta);
% y = r*sin(Omega*t)- L *cos(theta)
% plot(x,y, '.')
plot(theta,theta_point, '.r')
ms = 2
lw = 2
fs = 16
figure 
plot(theta,theta_point, '.r','MarkerSize',1.5)
grid on
xlabel('$\theta$ [rad]','FontSize',20,'Interpreter','latex')
ylabel(['$\dot \theta$ [rad/s]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)

%Conditions initiales
% tFin =85110.2012065732 % 100000 fois la periode 0.851102012065732
% r = 0.09
% Omega = 7.382411530116701
% kappa = 0.1
% m = 0.12
% g = 9.81
% L = 0.18
% theta0 = 0
% thetadot0 = 0.0
% nsteps = 100000000 %100000 periodes fois 1000 iteration par periode
% output = poincarre.out
% sampling = 1000 
