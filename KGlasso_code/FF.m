function [Frob_error_ML_inv_final, Frob_error_ML_cov_final, A_mat, B_mat] = FF(SCM,p,f,n,niter, tol,A0,B0,X0,Y0,reg)
%FF Function implements the flip-flop algorithm.
% 
% Last updated: Sept. 18, 2011
% 

% global SCM;
% tol = 1e-8;


% initial estimates of X, Y matrices
X_mat = eye(p);
Y_mat = eye(f);

% Likelihood = zeros(niter,1);
% Frob_error_ML_vec = zeros(niter,1);
% Frob_error_ML_cov_vec = zeros(niter,1);

% Theta_est = kron(X_mat,Y_mat); % initial guess

% temp1 = trace(Theta_est*SCM) - log(det(Theta_est));
% temp2 = norm(CovMat_inv-Theta_est,'fro')/norm(CovMat_inv,'fro');

% % compute ML solution (as in Werner's paper)
% Cov_inv_current = kron(X_mat,Y_mat);
% abs_diff = 10000;

for ni=1:niter,
%     Theta_prev = Theta_est;
    X_prev = X_mat;
    Y_prev = Y_mat;
        
    % compute temporary matrix T_X
    T_X = zeros(f,f);
    for i1=1:p,
        for i2=1:p,
        	T_X = T_X + X_mat(i2,i1)*SCM((i1-1)*f+1:i1*f,(i2-1)*f+1:i2*f);
        end
    end
    T_X = (T_X + reg*trace(X_mat)*eye(f))/p;
    B_hat = T_X;
    B_hat_inv = inv(T_X);    
    Y_mat = B_hat_inv;
    %threshold
    %Y_mat = arrayfun(@thresh, Y_mat);
    % compute T_Y
    T_Y = zeros(p,p);
    for j1=1:f,
        for j2=1:f,
        	T_Y = T_Y + Y_mat(j2,j1)*SCM(j1:f:end,j2:f:end);
        end
    end
    T_Y = (T_Y + reg*trace(Y_mat)*eye(p))/f;
    A_hat = T_Y;
    A_hat_inv = inv(T_Y);
    X_mat = A_hat_inv;
    %threshold
    %X_mat = arrayfun(@thresh, X_mat);
    
%     Theta_est = kron(X_mat, Y_mat);
    
%     Frob_diff = norm(Theta_prev-Theta_est,'fro');
    Frob_diff = sqrt(computeFrob(X_mat,Y_mat,X_prev,Y_prev));
    Frob_diff = Frob_diff./norm(X_prev,'fro')*norm(Y_prev,'fro');
    
%     Frob_diff, pause
    
    if Frob_diff<tol,
%         disp(['Exited ML solution at FF iteration = ' num2str(ni)]);
%         pause
        break;
    end

%     % store objective value at end of each iteration and compute error in
%     % overall covariance matrix (to plot later)
%     Theta_est = kron(X_mat,Y_mat);

%     lik_est = trace(Theta_est*SCM) - log(det(Theta_est));
%     Likelihood(ni) = lik_est;
%     Frob_error_ML_vec(ni) = norm(CovMat_inv-Theta_est,'fro')^2/norm(CovMat_inv,'fro')^2; % inverse
%     Frob_error_ML_cov_vec(ni) = norm(CovMat-kron(A0,B_hat),'fro')^2/norm(CovMat,'fro')^2; % forward
    
end

% lik_optimal = trace(CovMat_inv*SCM) - log(det(CovMat_inv));
% Likelihood = [temp1; Likelihood] - lik_optimal*ones(niter+1,1);
% Frob_error_ML_vec = [temp2; Frob_error_ML_vec];

% Theta_0 = kron(X0,Y0);
% Sigma_0 = kron(A0,B0);
% Frob_error_ML_inv_final = norm(Theta_0-Theta_est,'fro')^2/norm(Theta_0,'fro')^2; % inverse
% Frob_error_ML_cov_final = norm(Sigma_0-kron(A_hat,B_hat),'fro')^2/norm(Sigma_0,'fro')^2; % forward

% whos
    
temp1 = (norm(X0,'fro')*norm(Y0,'fro'))^2;
temp2 = (norm(A0,'fro')*norm(B0,'fro'))^2;

%thresholding
X_mat = arrayfun(@thresh, X_mat);
Y_mat = arrayfun(@thresh, Y_mat);

A_hat = inv(X_mat)
B_hat = inv(Y_mat)

Frob_error_ML_inv_final = computeFrob(X_mat,Y_mat,X0,Y0)/temp1; % inverse
Frob_error_ML_cov_final = computeFrob(A_hat,B_hat,A0,B0)/temp2; % forward


% Frob_error_ML_inv_final
% Frob_error_ML_cov_final
% pause



% % store output
% % ni_ML = ni;
% A_hat_ML = inv(X_mat);
% B_hat_ML = inv(Y_mat);
% CovMat_est_ML = kron(A_hat_ML, B_hat_ML);
% CovMat_est_inv_ML = Cov_inv_current;



A_mat = A_hat;
B_mat = B_hat;



end
