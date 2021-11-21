%caos r = 0.14 omega0=0/1e-10 tfin 100s   κ = 0.01, theta_point_0= 0
data = load('c_3_caos_theta=0.out');
t = data(:,1);
theta1= data(:,2);
data = load('c_3_caos_theta=1e-10.out');
theta2= data(:,2);
err= abs(theta1-theta2);



ms = 2
lw = 2
fs = 16
figure 
plot(t, theta1, '.r','LineWidth',lw,'MarkerSize',ms)
hold on
plot(t, theta2, '.r','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('$t$ [s]','FontSize',20,'Interpreter','latex')
ylabel(['$\theta$ [rad]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)
%%après un petit temps les 2 courbe ne suive pas la meme courbe


figure 
plot(t, err,'.r','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('$t$ [s]','FontSize',20,'Interpreter','latex')
ylabel(['Erreur sur $\theta$ [rad]]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)
set(gca,'YScale', 'log')
%erreur augment presque linéairement au début et après se stabilise



%no_caos r = 0.01 omega0=0/1e-10 tfin 100s
%caos r = 0.14 omega0=0/1e-10 tfin 100s
data = load('c_3_no_caos_theta=0.out');
t = data(:,1);
theta1= data(:,2);
data = load('c_3_no_caos_theta=1e-10.out');
theta2= data(:,2);
err= abs(theta1-theta2);
figure 
plot(t, theta1, '.r','LineWidth',lw,'MarkerSize',ms)
hold on
plot(t, theta2, '.r','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('$t$ [s]','FontSize',20,'Interpreter','latex')
ylabel(['$\theta$ [rad]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)
%% oscillation que se stabilisent et restent costantes %% exactement la même courbe pour theta 1 ou theta 2


figure 
plot(t, err,'.r','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('$t$ [s]','FontSize',20,'Interpreter','latex')
ylabel(['Erreur sur $\theta$ [rad]]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)
set(gca,'YScale', 'log') 
%% theta_1 - theta_2 très très petit
