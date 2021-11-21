lw = 2
fs = 16
ms = 11
theta_ini = [2.80]
E_ini = [0.199653145098326]
t = [2.5]
delta_t =[40 60 80 160 320] 
E_fin = [0.199881648336226 0.199754354672064 0.199709922794652 0.19966730203849 0.199656681967132] %0.206244518797253
figure 
plot((t./delta_t).^2, E_fin-E_ini, '+r','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel(['$\Delta t ^2$ [s$^2$]'],'FontSize',20,'Interpreter','latex')
ylabel(['$Erreur sur E_{mec}$ [J]'],'Fontsize', 20,'Interpreter','latex')
set(gca,'FontSize',fs)

% fit y = 0.05852*x - 2.002e-07 %% noter que en générale l'erreur sur
% l'E_mec est très petit
%E_initiale = [0.199653145098326]
% t_simulation = [2.5]
% nsteps =[40 60 80 160 320]
%r = 0 k = 0 
% Omega = 7.382411530116701