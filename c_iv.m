% r = 0.14 Omega = 7.382411530116701 theta_0 = 0 theta_point_0 = 0
%nstep = 10000 20000 40000 80000 160000
%%les trajectoire après seulement 1.5 seconds sont très different
%

r = 0.14
L = 0.18
Omega = 7.382411530116701
data = load('c_4_nsteps=10000.out');
theta = data(:,2);
t = data(:,1);
x = r*cos(Omega*t) + L*sin(theta);
y = r*sin(Omega*t)- L *cos(theta)
ms = 2
lw = 2
fs = 16
figure 
plot(x, y, 'm','MarkerSize',ms)
grid on
xlabel('$t$ [s]','FontSize',20,'Interpreter','latex')
ylabel(['$\theta$ [rad]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)
hold on
data = load('c_4_nsteps=20000.out');
theta = data(:,2);
t = data(:,1);
x = r*cos(Omega*t) + L*sin(theta);
y = r*sin(Omega*t)- L *cos(theta)
plot(x, y, 'b','MarkerSize',ms)
% plot(t, theta, '.b','MarkerSize',ms)
% hold on
% % data = load('c_4_nsteps=40000.out');
% % theta = data(:,2);
% % t = data(:,1);
% % x = r*cos(Omega*t) + L*sin(theta);
% % y = r*sin(Omega*t)- L *cos(theta)
% % plot(x, y, '.g','MarkerSize',ms)

% plot(t, theta, '.r','MarkerSize',ms)
% hold on
% data = load('c_4_nsteps=80000.out');
% theta = data(:,2);
% t = data(:,1);
% plot(t, theta, '.g','MarkerSize',ms)
% 
% hold on
% data = load('c_4_nsteps=160000.out');
% theta = data(:,2);
% t = data(:,1);
% plot(t, theta, '.b','MarkerSize',ms)
% 
