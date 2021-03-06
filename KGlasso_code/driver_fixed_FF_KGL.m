% Script file compares performance between: iterative KGL, 
% FF, naive Glasso.
% 
% Last updated: May 3, 2012
%
clear all, close all;
randn('state',1);
rand('state',22);

boolGlasso = 0;
boolSCM = 0;
boolLarge = 1;

% parameters
niterKGL   = 20;
niterFF    = 20;
ntrials    = 5; % number of Monte Carlo runs
ts = datestr(now);
timestamp = [ts(1:11) '-' ts(13:14) '-' ts(16:17)];

p=10; f=20; % ex. 1,2
%p = 30; f = 60;
% p=100; f=100;

tol = 1e-5; % convergence tolerance for each iterative method (FF,KGL)

% randomly choose A0, B0
[A0,X0] = genSparseER(p);
[B0,Y0] = genSparseER(f);

% display matrices
figure;
subplot(121), imagesc(X0), axis square, colormap(gray), colormap(1-colormap), title('X_0');
subplot(122), imagesc(Y0), axis square, colormap(gray), colormap(1-colormap), title('Y_0');
% subplot(133), imagesc(kron(X0,Y0)), axis square, colormap(gray), title('\Theta_0=X_0 \otimes Y_0');
files.X0 = X0;
files.Y0 = Y0;
files.A0 = A0;
files.B0 = B0;
saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_matrices.fig']);
X_off_diag = X0-diag(diag(X0));
s_X0 = sum(X_off_diag(:)~=0);
Y_off_diag = Y0-diag(diag(Y0));
s_Y0 = sum(Y_off_diag(:)~=0);
sparsity_string  = sprintf('s_X0 = %d, s_Y0 = %d', s_X0, s_Y0);
sparsity_stringx = sprintf('percentage of zeros on off-diagonal of X0 = %0.2f', 100*(1-s_X0/(p^2-p)));
sparsity_stringy = sprintf('percentage of zeros on off-diagonal of Y0 = %0.2f', 100*(1-s_Y0/(f^2-f)));
disp(sparsity_string);
disp(sparsity_stringx);
disp(sparsity_stringy);
pause % ok good

len_n = 10;
% n_vec = floor(logspace(1,5,len_n)).' % large-sample run
% n_vec = floor(logspace(0,2,len_n)).'
n_vec = 3*(1:len_n)
reg = 2*(1 + log(n_vec)).^(-1)
% disp(n_vec)

zero_vec = zeros(length(n_vec),1);

Frob_error_inv_KGL_vec = zero_vec;
Frob_error_inv_FFP_vec = zero_vec;
Frob_error_inv_FF_vec = zero_vec;
Frob_error_inv_regFF_vec = zero_vec;
Frob_error_inv_Glasso_vec = zero_vec;
Frob_error_inv_SCM_vec = zero_vec;

Frob_error_inv_KGL_std_vec = zero_vec;
Frob_error_inv_FF_std_vec = zero_vec;
Frob_error_inv_regFF_std_vec = zero_vec;
Frob_error_inv_Glasso_std_vec = zero_vec;
Frob_error_inv_SCM_std_vec = zero_vec;

Frob_error_cov_KGL_vec = zero_vec;
Frob_error_cov_FF_vec = zero_vec;
Frob_error_cov_regFF_vec = zero_vec;
Frob_error_cov_Glasso_vec = zero_vec;
Frob_error_cov_SCM_vec = zero_vec;

Frob_error_cov_KGL_std_vec = zero_vec;
Frob_error_cov_FF_std_vec = zero_vec;
Frob_error_cov_regFF_std_vec = zero_vec;
Frob_error_cov_Glasso_std_vec = zero_vec;
Frob_error_cov_SCM_std_vec = zero_vec;

count = 0;

% CovMat = kron(A0,B0);
% CovMat_inv = kron(X_true,Y_true);

hh = waitbar(0,'Please wait...');
for kk=1:length(n_vec),
    temp = zeros(ntrials,6);
    tempFFP = zeros(ntrials,2);
    tempGlasso = zeros(ntrials,2);
    tempSCM = zeros(ntrials,2);
    for nn=1:ntrials,
        % generate Gaussian data with Kronecker structure
        [SCM] = dataGenFixed_v2(n_vec(kk),A0,B0);
                
        % KGL - Kronecker graphical lasso
        cx = 0.4;
        cy = 0.4;
        %cx = 1
        %cy = 1
        
        tic;
        [err_inv_KGL, err_cov_KGL] = KGL_iterative(SCM,p,f,n_vec(kk),cx,cy,A0,B0,X0,Y0,niterKGL,tol);
        t=toc;
        disp(['KGL timing = ' num2str(t) ' secs']);
        
        
        % FF - standard flip-flop
        
        if (n_vec(kk)>3),
            tic;
            [err_inv_FF, err_cov_FF] = FF(SCM,p,f,n_vec(kk),niterFF,tol,A0,B0,X0,Y0,0);
            t=toc;
            disp(['FF timing = ' num2str(t) ' secs']);
        else
            err_inv_FF = 0
            err_cov_FF = 0
        end
        
        tic;
        [err_inv_regFF, err_cov_regFF] = FF(SCM,p,f,n_vec(kk),niterFF,tol,A0,B0,X0,Y0,reg(kk));
        t=toc;
        disp(['regFF timing = ' num2str(t) ' secs']);
        
        % store results in temporary vector
        temp_vec = [err_inv_KGL, err_inv_FF, ...
                    err_cov_KGL, err_cov_FF, err_inv_regFF, err_cov_regFF];
        temp(nn,:) = temp_vec;
        
        % display error during run...
        [sqrt(mean(temp(nn,1))) sqrt(mean(temp(nn,2))) sqrt(mean(temp(nn,3))) sqrt(mean(temp(nn,4))) sqrt(mean(temp(nn,5))) sqrt(mean(temp(nn,6)))]

        % Standard GLASSO
        if boolGlasso,
            c = 1.1;
            
            lambdaGlasso = c*sqrt(log(p*f)/n_vec(kk));
            disp('Entering Naive Glasso...');
            tic;
            [err_inv_Glasso, err_cov_Glasso] = EstGlasso(SCM,lambdaGlasso,A0,B0,X0,Y0);
            t=toc;
            disp(['Naive Glasso timing = ' num2str(t) ' secs']);
            disp('Exited Naive Glasso!');
            tempGlasso(nn,:) = [err_inv_Glasso err_cov_Glasso]
        end
        
        % Standard SCM
        if boolSCM,
            err_inv_SCM = norm(pinv(SCM)-kron(X0,Y0),'fro')^2/(norm(kron(X0,Y0),'fro')^2);
            err_cov_SCM = norm(SCM-kron(A0,B0),'fro')^2/(norm(kron(A0,B0),'fro')^2);
            tempSCM(nn,:) = [err_inv_SCM err_cov_SCM];
        end
        
        count = count + 1;
        waitbar(count/(length(n_vec)*ntrials),hh,sprintf('%3.2f%% Done', count/(length(n_vec)*ntrials)*100));
    end

    % normalized metrics
    Frob_error_inv_KGL_vec(kk)     = sqrt(mean(temp(:,1)));
    Frob_error_inv_FF_vec(kk)      = sqrt(mean(temp(:,2)));
    Frob_error_inv_regFF_vec(kk)      = sqrt(mean(temp(:,5)));
    Frob_error_cov_KGL_vec(kk)     = sqrt(mean(temp(:,3)));
    Frob_error_cov_FF_vec(kk)      = sqrt(mean(temp(:,4)));
    Frob_error_cov_regFF_vec(kk)      = sqrt(mean(temp(:,6)));
    % compute std
    Frob_error_inv_KGL_std_vec(kk) = std(sqrt(temp(:,1)));
    Frob_error_inv_FF_std_vec(kk)  = std(sqrt(temp(:,2)));
    Frob_error_cov_KGL_std_vec(kk) = std(sqrt(temp(:,3)));
    Frob_error_cov_FF_std_vec(kk)  = std(sqrt(temp(:,4)));
    Frob_error_inv_regFF_std_vec(kk)  = std(sqrt(temp(:,5)));
    Frob_error_cov_regFF_std_vec(kk)  = std(sqrt(temp(:,6)));
    
    if boolSCM,
        Frob_error_inv_SCM_vec(kk) = sqrt(mean(tempSCM(:,1)));
        Frob_error_cov_SCM_vec(kk) = sqrt(mean(tempSCM(:,2)));
    end
    
    if boolGlasso,
        Frob_error_inv_Glasso_vec(kk) = sqrt(mean(tempGlasso(:,1)));
        Frob_error_cov_Glasso_vec(kk) = sqrt(mean(tempGlasso(:,2)));

        Frob_error_inv_Glasso_std_vec(kk) = std(sqrt(tempGlasso(:,1)));
        Frob_error_cov_Glasso_std_vec(kk) = std(sqrt(tempGlasso(:,2)));
    end 
    
%     inv_cov_err = [Frob_error_inv_KGL_vec Frob_error_inv_FF_vec Frob_error_inv_Glasso_vec].'
%     cov_err = [Frob_error_cov_KGL_vec Frob_error_cov_FF_vec Frob_error_cov_Glasso_vec].'   
%     pause
end
close(hh);

disp(reg)

%% PLOTS

figure;
loglog(n_vec, Frob_error_inv_FF_vec,'-ro'), hold on, loglog(n_vec, Frob_error_inv_regFF_vec,'-gd'), hold on, loglog(n_vec, Frob_error_inv_KGL_vec,'-cx'),...
    xlabel('Sample size (n)'), ylabel('Normalized Frobenius error in inverse cov. matrix'), ylim([1e-2 1e1]), legend('FF','regFF','KGL'), axis square, ...
    title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_inv_all.fig']);

% figure;
% loglog(n_vec, Frob_error_cov_FF_vec,'-ro'), hold on, loglog(n_vec, Frob_error_cov_regFF_vec,'-go'),hold on, loglog(n_vec, Frob_error_cov_KGL_vec,'-cx'),...
%     xlabel('Sample size (n)'), ylabel('Normalized Frobenius error in cov. matrix'), ylim([1e-2 1e1]), legend('FF','regFF','KGL'), axis square, ...
%     title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
% saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_cov_all.fig']);

% Standard Dev. plots
figure;
errorbar(n_vec, Frob_error_inv_FF_vec,2*Frob_error_inv_FF_std_vec,'Color','r','Marker','o','LineStyle','-'), ...
    hold on, errorbar(n_vec, Frob_error_inv_regFF_vec,2*Frob_error_inv_regFF_std_vec,'Color','g','Marker','o','LineStyle','-'), ...
    hold on, errorbar(n_vec, Frob_error_inv_KGL_vec,2*Frob_error_inv_KGL_std_vec,'Color','c','Marker','x','LineStyle','-'),...
    set(gca,'xscale','log'), set(gca,'yscale','log'),...
    xlabel('Sample size (n)'), ylabel('Normalized Frobenius error in inverse cov. matrix'), ylim([1e-2 1e1]), legend('FF','regFF','KGL'), axis square, ...
    title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_inv_std_all.fig']);

% figure;
% errorbar(n_vec, Frob_error_cov_FF_vec,2*Frob_error_cov_FF_std_vec,'Color','r','Marker','o','LineStyle','-'), ...
%     hold on, errorbar(n_vec, Frob_error_cov_regFF_vec,2*Frob_error_cov_regFF_std_vec,'Color','g','Marker','o','LineStyle','-'), ...
%     hold on, errorbar(n_vec, Frob_error_cov_KGL_vec,2*Frob_error_cov_KGL_std_vec,'Color','c','Marker','x','LineStyle','-'),...
%     set(gca,'xscale','log'), set(gca,'yscale','log'),...
%     xlabel('Sample size (n)'), ylabel('Normalized Frobenius error in cov. matrix'), ylim([1e-2 1e1]), legend('FF','regFF','KGL'), axis square, ...
%     title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
% saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_cov_std_all.fig']);

if boolGlasso,
    figure;
    loglog(n_vec, Frob_error_inv_KGL_vec,'-bv'), hold on, loglog(n_vec, Frob_error_inv_FF_vec,'-ro'), hold on, loglog(n_vec, Frob_error_inv_Glasso_vec,'-ks'),...
        xlabel('Sample size (n)'), ylabel('Normalized Frob. error in precision matrix'), ylim([1e-2 1e1]), legend('KGL','FF','Glasso'), axis square, ...
        title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
    saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_inv_withGlasso.fig']);

    figure;
    loglog(n_vec, Frob_error_cov_KGL_vec,'-bv'), hold on, loglog(n_vec, Frob_error_cov_FF_vec,'-ro'), hold on, loglog(n_vec, Frob_error_cov_Glasso_vec,'-ks'),...
        xlabel('Sample size (n)'), ylabel('Normalized Frob. error in covariance matrix'), ylim([1e-2 1e1]), legend('KGL','FF','Glasso'), axis square, ...
        title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
    saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_cov_withGlasso.fig']);
    
    if boolSCM,
        figure;
        loglog(n_vec, Frob_error_inv_KGL_vec,'-bv'), hold on, loglog(n_vec, Frob_error_inv_FF_vec,'-ro'), ...
            hold on, loglog(n_vec, Frob_error_inv_Glasso_vec,'-ks'), hold on, loglog(n_vec, Frob_error_inv_SCM_vec,'-gd'),...
            xlabel('Sample size (n)'), ylabel('Normalized Frob. error in precision matrix'), ylim([1e-2 1e1]), legend('KGL','FF','Glasso','SCM'), axis square, ...
            title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
        saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_inv_withGlasso_withSCM.fig']);

        figure;
        loglog(n_vec, Frob_error_cov_KGL_vec,'-bv'), hold on, loglog(n_vec, Frob_error_cov_FF_vec,'-ro'), ...
            hold on, loglog(n_vec, Frob_error_cov_Glasso_vec,'-ks'), hold on, loglog(n_vec, Frob_error_cov_SCM_vec,'-gd'),...
            xlabel('Sample size (n)'), ylabel('Normalized Frob. error in covariance matrix'), ylim([1e-2 1e1]), legend('KGL','FF','Glasso','SCM'), axis square, ...
            title(['Max number of iter. = ' num2str(niterKGL) ', Trials = ' num2str(ntrials) ', (p,f)=(' num2str(p) ',' num2str(f) ')']);
        saveas(gcf, ['./SimSynthetic3/KGL_FF_sim_' timestamp '_Frob_cov_withGlasso_withSCM.fig']);
    end
end


% save matrices in mat file
save(['./SimSynthetic3/KGL_FF_sim_' timestamp '.mat'], '-struct', 'files');

% store all variable for later processing
save(['./SimSynthetic3/KGL_FF_sim_' timestamp 'all.mat']);

