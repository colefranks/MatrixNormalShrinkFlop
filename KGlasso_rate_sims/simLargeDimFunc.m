function [err_empirical err_predicted] = simLargeDimFunc(alpha, n_vec, scaling)
%SIMLARGEDIMFUNC Functions compares high dimensional performance of FF and
%KGLASSO.

% specs
niterKGL   = 20;
niterFF    = 20;
boolLarge  = 1;
tol        = 1e-3; % convergence tolerance for each iterative method
ntrials    = 20;

p_vec = ceil(scaling*n_vec.^alpha);
[n_vec(:) p_vec(:)]
% KGL_pred = sqrt(p_vec.*log(max(p_vec,n_vec))./n_vec);
% figure, plot(n_vec, KGL_pred)
% pause

zero_vec = zeros(length(n_vec),1);
Frob_error_inv_KGL_vec = zero_vec;
Frob_error_inv_KGL_std_vec = zero_vec;
Frob_error_cov_KGL_vec = zero_vec;
Frob_error_cov_KGL_std_vec = zero_vec;
Frob_error_inv_FF_vec = zero_vec;
Frob_error_inv_FF_std_vec = zero_vec;
Frob_error_cov_FF_vec = zero_vec;
Frob_error_cov_FF_std_vec = zero_vec;

count = 0;
hh = waitbar(0,'Please wait...');
for i=1:length(n_vec),
    p = p_vec(i);
    f = p;
    
    A0 = eye(p);
    B0 = eye(f);
    X0 = eye(p);
    Y0 = eye(f);
%     [trash,trash,A0,B0,X0,Y0,Sigma_0,Theta_0] = dataGen(p,f,2,0);
%     [n_vec(i), p_vec(i)]
%     
%     X_off_diag = X0-diag(diag(X0));
%     s_X0 = sum(X_off_diag(:)~=0);
%     Y_off_diag = Y0-diag(diag(Y0));
%     s_Y0 = sum(Y_off_diag(:)~=0);
%     sparsity_string  = sprintf('s_X0 = %d, s_Y0 = %d', s_X0, s_Y0);
%     sparsity_stringx = sprintf('percentage of zeros on off-diagonal of X0 = %0.2f', 100*(1-s_X0/(p^2-p)));
%     sparsity_stringy = sprintf('percentage of zeros on off-diagonal of Y0 = %0.2f', 100*(1-s_Y0/(f^2-f)));
%     disp(sparsity_string);
%     disp(sparsity_stringx);
%     disp(sparsity_stringy);
% %     pause % ok good

    temp = zeros(ntrials,4);
    for nn=1:ntrials,
        [SCM] = dataGenFixed_v2(n_vec(i),A0,B0);
        
        tic;
        [err_inv_FF, err_cov_FF] = FF(SCM,p,f,n_vec(i),niterFF,tol,A0,B0,X0,Y0);
        err_inv_FF = err_inv_FF*(norm(X0,'fro')*norm(Y0,'fro'))^2; % de-normalize MSE (dimension is varying here)
        err_cov_FF = err_cov_FF*(norm(A0,'fro')*norm(B0,'fro'))^2;
        t=toc;
        disp(['FF timing = ' num2str(t) ' secs']);

        
        cx = 0.4;
        cy = 0.4;
        tic;
        [err_inv_KGL, err_cov_KGL] = KGL_iterative(SCM,p,f,n_vec(i),cx,cy,A0,B0,X0,Y0,niterKGL,tol);
        err_inv_KGL = err_inv_KGL*(norm(X0,'fro')*norm(Y0,'fro'))^2; % de-normalize MSE (dimension is varying here)
        err_cov_KGL = err_cov_KGL*(norm(A0,'fro')*norm(B0,'fro'))^2;
        t=toc;
        disp(['KGL timing = ' num2str(t) ' secs']);
        
        
        temp_vec = [err_inv_KGL err_cov_KGL err_inv_FF err_cov_FF];
        temp(nn,:) = temp_vec;    
        
        % display during run...
        [sqrt(mean(temp(nn,2))) sqrt(mean(temp(nn,4)))]
        
        count = count + 1;
        waitbar(count/(length(n_vec)*ntrials),hh,sprintf('%3.2f%% Done', count/(length(n_vec)*ntrials)*100));
    end
    
    % normalized error metrics
    Frob_error_inv_KGL_vec(i)     = sqrt(mean(temp(:,1)));
    Frob_error_cov_KGL_vec(i)     = sqrt(mean(temp(:,2)));
    Frob_error_inv_KGL_std_vec(i) = std(sqrt(temp(:,1)));
    Frob_error_cov_KGL_std_vec(i) = std(sqrt(temp(:,2)));
    Frob_error_inv_FF_vec(i)      = sqrt(mean(temp(:,3)));
    Frob_error_cov_FF_vec(i)      = sqrt(mean(temp(:,4)));
    Frob_error_inv_FF_std_vec(i)  = std(sqrt(temp(:,3)));
    Frob_error_cov_FF_std_vec(i)  = std(sqrt(temp(:,4)));    
end
close(hh);


%% plot results

% covariance matrix
KGL_pred = sqrt(p_vec.*log(max(p_vec,n_vec))./n_vec);
FF_pred  = sqrt(p_vec.^2./n_vec);
KGL_pred = KGL_pred./KGL_pred(end)*Frob_error_cov_KGL_vec(end);
FF_pred  = FF_pred./FF_pred(end)*Frob_error_cov_FF_vec(end);
err_predicted.cov.FF = FF_pred;
err_predicted.cov.KGL = KGL_pred;
err_empirical.cov.FF = Frob_error_cov_FF_vec;
err_empirical.cov.KGL = Frob_error_cov_KGL_vec;

% precision matrix
KGL_pred = sqrt(p_vec.*log(max(p_vec,n_vec))./n_vec);
FF_pred = sqrt(p_vec.^2./n_vec);
KGL_pred = KGL_pred./KGL_pred(end)*Frob_error_inv_KGL_vec(end);
FF_pred  = FF_pred./FF_pred(end)*Frob_error_inv_FF_vec(end);
err_predicted.inv.FF = FF_pred;
err_predicted.inv.KGL = KGL_pred;
err_empirical.inv.FF = Frob_error_inv_FF_vec;
err_empirical.inv.KGL = Frob_error_inv_KGL_vec;




end

