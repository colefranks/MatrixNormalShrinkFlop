% Script simulates a family of curves.
close all, clear all;

ts = datestr(now);
timestamp = [ts(1:11) '-' ts(13:14) '-' ts(16:17)];

alpha = [0.1 0.2 0.3];
n_vec = floor(linspace(5,1000,10));
% n_vec = 1000
constant = 8;
ceil(constant*n_vec.^alpha(3))
pause

[err_frob_empirical_1 err_frob_predicted_1] = simLargeDimFunc(alpha(1),n_vec,constant);
[err_frob_empirical_2 err_frob_predicted_2] = simLargeDimFunc(alpha(2),n_vec,constant);
[err_frob_empirical_3 err_frob_predicted_3] = simLargeDimFunc(alpha(3),n_vec,constant);

% save files in MAT file
save(['./SimNew/KGlasso_FF_sim_' timestamp '_all.mat']);


%% Plots results
figure(1); % FF empirical
plot(n_vec,err_frob_empirical_1.inv.FF,'-bo');
hold on, plot(n_vec,err_frob_empirical_2.inv.FF,'-rx');
hold on, plot(n_vec,err_frob_empirical_3.inv.FF,'-g+');
legend(['\alpha = ' num2str(alpha(1))],['\alpha = ' num2str(alpha(2))],['\alpha = ' num2str(alpha(3))]);
title('FF'), ylabel('Empirical MSE for Precision Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_inv_FF.fig']);

figure(2); % KGL empirical
plot(n_vec,err_frob_empirical_1.inv.KGL,'-bo');
hold on, plot(n_vec,err_frob_empirical_2.inv.KGL,'-rx');
hold on, plot(n_vec,err_frob_empirical_3.inv.KGL,'-g+');
legend(['\alpha = ' num2str(alpha(1))],['\alpha = ' num2str(alpha(2))],['\alpha = ' num2str(alpha(3))]);
title('KGlasso'), ylabel('Empirical MSE for Precision Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_inv_KGL.fig']);

figure(3); % FF empirical + prediction
semilogy(n_vec,err_frob_empirical_1.inv.FF,'-bo'), hold on, semilogy(n_vec,err_frob_predicted_1.inv.FF,'-bd');
hold on, semilogy(n_vec,err_frob_empirical_2.inv.FF,'-rx'), hold on, semilogy(n_vec,err_frob_predicted_2.inv.FF,'-rs');
hold on, semilogy(n_vec,err_frob_empirical_3.inv.FF,'-g+'), hold on, semilogy(n_vec,err_frob_predicted_3.inv.FF,'-gh');
legend(['Emp. \alpha = ' num2str(alpha(1))],['Pred. \alpha = ' num2str(alpha(1))],['Emp. \alpha = ' num2str(alpha(2))],['Pred. \alpha = ' num2str(alpha(2))],['Emp. \alpha = ' num2str(alpha(3))],['Pred. \alpha = ' num2str(alpha(3))]);
title('FF'), ylabel('Empirical MSE for Precision Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_inv_FF_pred.fig']);

figure(4); % KGL empirical + prediction
semilogy(n_vec,err_frob_empirical_1.inv.KGL,'-bo'), hold on, semilogy(n_vec,err_frob_predicted_1.inv.KGL,'-bd');
hold on, semilogy(n_vec,err_frob_empirical_2.inv.KGL,'-rx'), hold on, semilogy(n_vec,err_frob_predicted_2.inv.KGL,'-rs');
hold on, semilogy(n_vec,err_frob_empirical_3.inv.KGL,'-g+'), hold on, semilogy(n_vec,err_frob_predicted_3.inv.KGL,'-gh');
legend(['Emp. \alpha = ' num2str(alpha(1))],['Pred. \alpha = ' num2str(alpha(1))],['Emp. \alpha = ' num2str(alpha(2))],['Pred. \alpha = ' num2str(alpha(2))],['Emp. \alpha = ' num2str(alpha(3))],['Pred. \alpha = ' num2str(alpha(3))]);
title('KGlasso'), ylabel('Empirical MSE for Precision Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_inv_KGL_pred.fig']);

figure(5); % FF empirical + prediction
semilogy(n_vec,err_frob_empirical_1.cov.FF,'-bo'), hold on, semilogy(n_vec,err_frob_predicted_1.cov.FF,'-bd');
hold on, semilogy(n_vec,err_frob_empirical_2.cov.FF,'-rx'), hold on, semilogy(n_vec,err_frob_predicted_2.cov.FF,'-rs');
hold on, semilogy(n_vec,err_frob_empirical_3.cov.FF,'-g+'), hold on, semilogy(n_vec,err_frob_predicted_3.cov.FF,'-gh');
legend(['Emp. \alpha = ' num2str(alpha(1))],['Pred. \alpha = ' num2str(alpha(1))],['Emp. \alpha = ' num2str(alpha(2))],['Pred. \alpha = ' num2str(alpha(2))],['Emp. \alpha = ' num2str(alpha(3))],['Pred. \alpha = ' num2str(alpha(3))]);
title('FF'), ylabel('Empirical MSE for Covariance Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_cov_FF_pred.fig']);

figure(6); % KGL empirical + prediction
semilogy(n_vec,err_frob_empirical_1.cov.KGL,'-bo'), hold on, semilogy(n_vec,err_frob_predicted_1.cov.KGL,'-bd');
hold on, semilogy(n_vec,err_frob_empirical_2.cov.KGL,'-rx'), hold on, semilogy(n_vec,err_frob_predicted_2.cov.KGL,'-rs');
hold on, semilogy(n_vec,err_frob_empirical_3.cov.KGL,'-g+'), hold on, semilogy(n_vec,err_frob_predicted_3.cov.KGL,'-gh');
legend(['Emp. \alpha = ' num2str(alpha(1))],['Pred. \alpha = ' num2str(alpha(1))],['Emp. \alpha = ' num2str(alpha(2))],['Pred. \alpha = ' num2str(alpha(2))],['Emp. \alpha = ' num2str(alpha(3))],['Pred. \alpha = ' num2str(alpha(3))]);
title('KGlasso'), ylabel('Empirical MSE for Covariance Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_cov_KGL_pred.fig']);

figure(7); % summary of fig. 3,4
semilogy(n_vec,err_frob_empirical_1.inv.FF,':bo'), hold on, semilogy(n_vec,err_frob_predicted_1.inv.FF,':bd');
hold on, semilogy(n_vec,err_frob_empirical_2.inv.FF,':rx'), hold on, semilogy(n_vec,err_frob_predicted_2.inv.FF,':rs');
hold on, semilogy(n_vec,err_frob_empirical_3.inv.FF,':g+'), hold on, semilogy(n_vec,err_frob_predicted_3.inv.FF,':gh');
hold on, semilogy(n_vec,err_frob_empirical_1.inv.KGL,'-bo'), hold on, semilogy(n_vec,err_frob_predicted_1.inv.KGL,'-bd');
hold on, semilogy(n_vec,err_frob_empirical_2.inv.KGL,'-rx'), hold on, semilogy(n_vec,err_frob_predicted_2.inv.KGL,'-rs');
hold on, semilogy(n_vec,err_frob_empirical_3.inv.KGL,'-g+'), hold on, semilogy(n_vec,err_frob_predicted_3.inv.KGL,'-gh');
legend(['FF: Emp. \alpha = ' num2str(alpha(1))],['FF: Pred. \alpha = ' num2str(alpha(1))],['FF: Emp. \alpha = ' num2str(alpha(2))],['FF: Pred. \alpha = ' num2str(alpha(2))],['FF: Emp. \alpha = ' num2str(alpha(3))],['FF: Pred. \alpha = ' num2str(alpha(3))],['KGL: Emp. \alpha = ' num2str(alpha(1))],['KGL: Pred. \alpha = ' num2str(alpha(1))],['KGL: Emp. \alpha = ' num2str(alpha(2))],['KGL: Pred. \alpha = ' num2str(alpha(2))],['KGL: Emp. \alpha = ' num2str(alpha(3))],['KGL: Pred. \alpha = ' num2str(alpha(3))]);
title('KGlasso & FF'), ylabel('Empirical & Predicted MSE for Precision Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_inv_KGL_FF_pred.fig']);







