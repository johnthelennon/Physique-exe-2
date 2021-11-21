% COndition initiales tf in = 60s, avec r = 0.005m, κ = 0.01,theta0 = 0, theta_point0 = 0,


Omega = [6.5 6.7 6.9 7.1 7.2 7.3 7.382411530116701 7.5 7.6 7.8 8.0 8.2 6.8 6.85 6.88]  %6.865 6.86];

data = load('c_i_Omega=6.5.out');
theta = data(:,2);
M(1)= max(abs(theta));


data = load('c_i_Omega=6.7.out');
theta = data(:,2);
M(2)= max(abs(theta));


data = load('c_i_Omega=6.9.out');
theta = data(:,2);
M(3)= max(abs(theta));


data = load('c_i_Omega=7.1.out');
theta = data(:,2);
M(4)= max(abs(theta));


data = load('c_i_Omega=7.2.out');
theta = data(:,2);
M(5)= max(abs(theta));


data = load('c_i_Omega=7.3.out');
theta = data(:,2);
M(6)= max(abs(theta));


data = load('c_i_Omega=7.382.out');
theta = data(:,2);
M(7)= max(abs(theta));


data = load('c_i_Omega=7.5.out');
theta = data(:,2);
M(8)= max(abs(theta));


data = load('c_i_Omega=7.6.out');
theta = data(:,2);
M(9)= max(abs(theta));


data = load('c_i_Omega=7.8.out');
theta = data(:,2);
M(10)= max(abs(theta));


data = load('c_i_Omega=8.0.out');
theta = data(:,2);
M(11)= max(abs(theta));

data = load('c_i_Omega=8.2.out');
theta = data(:,2);
M(12)= max(abs(theta));

data = load('c_i_Omega=6.8.out');
theta = data(:,2);
M(13)= max(abs(theta));

data = load('c_i_Omega=6.85.out');
theta = data(:,2);
M(14)= max(abs(theta));

data = load('c_i_Omega=6.88.out');
theta = data(:,2);
M(15)= max(abs(theta));
% 
% data = load('c_i_Omega=6.865.out');
% theta = data(:,2);
% M(16)= max(abs(theta));
% 
% data = load('c_i_Omega=6.86.out');
% theta = data(:,2);
% M(17)= max(abs(theta));


ms = 11
lw = 2
fs = 16
figure 
plot(Omega, M,'+r','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel(['$\Omega$ [rad/s]'],'FontSize',20,'Interpreter','latex')
ylabel(['Max $\theta$ [J]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)