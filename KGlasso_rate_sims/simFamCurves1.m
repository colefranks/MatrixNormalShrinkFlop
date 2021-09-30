% Script simulates a family of curves.
close all, clear all;

ts = datestr(now);
timestamp = [ts(1:11) '-' ts(13:14) '-' ts(16:17)];

alpha = [0.6];
n_vec = floor(linspace(5,1000,10));
% n_vec = 1000
constant = 1;
ceil(constant*n_vec.^alpha)
pause

[err_frob_empirical err_frob_predicted] = simLargeDimFunc(alpha,n_vec,constant);

% save files in MAT file
save(['./SimNew/KGlasso_FF_sim_' timestamp '_all.mat']);


%% Plots results

figure(1);
semilogy(n_vec,err_frob_empirical.inv.FF,':bo'), hold on, semilogy(n_vec,err_frob_predicted.inv.FF,':bd');
hold on, semilogy(n_vec,err_frob_empirical.inv.KGL,'-rx'), hold on, semilogy(n_vec,err_frob_predicted.inv.KGL,'-rs');
legend(['FF: Emp. \alpha = ' num2str(alpha)],['FF: Pred. \alpha = ' num2str(alpha)],['KGL: Emp. \alpha = ' num2str(alpha)],['KGL: Pred. \alpha = ' num2str(alpha)]);
title('KGlasso & FF'), ylabel('Empirical & Predicted MSE for Precision Matrix'), xlabel('n');
saveas(gcf, ['./SimNew/KGlasso_FF_sim_' timestamp '_inv_KGL_FF_pred.fig']);
