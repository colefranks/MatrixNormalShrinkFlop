function [X_mat, Y_mat] = KGL_R(data,p,f,n,cx,cy,N_iter,tol)

%how do I want to pass this? Probably as a list of matrices.


%generate SCM from the data
SCM = zeros(p*f,p*f);
for m=1:n,
    SCM = SCM + data(:,m)*data(:,m)';
end
SCM = SCM/n;
                
% KGL - Kronecker graphical lasso
%cx = 0.4;
%cy = 0.4;


        
tic;
X_mat = eye(p); % starting point
Y_mat = eye(f);

% X_mat = genpd(p); % starting point
% Y_mat = eye(f);

% sparsify each factor
M = max(max(p,f),n);
lambda_initial = cy*sqrt(log(M)/(n*p));
lambda_final   = lambda_initial + cx*sqrt(log(M)/(n*f));
% lambda_x = cx*(1/sqrt(p)+1/sqrt(f))*sqrt(log(M)/n);
% lambda_y = cy*(1/sqrt(p)+1/sqrt(f))*sqrt(log(M)/n);

% progress = [];
for n1=1:N_iter,
%     Theta_prev = Theta_est;
    X_prev = X_mat;
    Y_prev = Y_mat;
    
%     disp(['Entering outer iteration ' num2str(n1)]);
    
%     % FF(start with A) + Sparsify B + FF(start with B) + sparsify A
%     [trash, trash, T_Y, T_X] = FF_init(SCM,p,f,n,N_iter_inner_max,tol_inner,A0,B0,X0,Y0,1,X_mat);

    % project onto B
    T_X = zeros(f,f);
    for i1=1:p,
        for i2=1:p,
            T_X = T_X + X_mat(i2,i1)*SCM((i1-1)*f+1:i1*f,(i2-1)*f+1:i2*f);
        end
    end
    T_X = T_X/p;
    
    if n1==1,
        lambda = lambda_initial;
    else
        lambda = lambda_final;
    end
%     lambda = lambda_y;

    % sparsify solution
    if lambda~=0,
%         tic;
        [B_hat_inv B_hat newtonIter] = quicGlasso(T_X, lambda, lambda, 1e-3, N_iter);
%         t=toc
%         newtonIter
%         B_hat_inv = glasso2(T_X, lambda);
%         B_hat = inv(B_hat_inv);
    else
        B_hat = T_X;
        B_hat_inv = inv(T_X);
    end
    inv_Y_mat = B_hat;
    Y_mat = B_hat_inv;
    
%     [trash, trash, T_Y, T_X] = FF_init(SCM,p,f,n,N_iter_inner_max,tol_inner,A0,B0,X0,Y0,0,Y_mat);
    
    % project onto A
    T_Y = zeros(p,p);
    for j1=1:f,
        for j2=1:f,
            T_Y = T_Y + Y_mat(j2,j1)*SCM(j1:f:end,j2:f:end);
        end
    end
    T_Y = T_Y/f;
    
%     lambda = lambda_x;
    
    % sparsify solution
    if lambda~=0,
%         tic
        [A_hat_inv A_hat newtonIter] = quicGlasso(T_Y, lambda, lambda, 1e-3, N_iter);
%         t=toc
%         newtonIter
%         A_hat_inv = glasso2(T_Y, lambda);
%         A_hat = inv(A_hat_inv);
    else
        A_hat = T_Y;
        A_hat_inv = inv(T_Y);
    end
    inv_X_mat = A_hat;
    X_mat = A_hat_inv;
    
    
%     % current estimate
%     Theta_est = kron(X_mat,Y_mat);
%     

    % check exit condition
    Frob_diff = sqrt(computeFrob(X_mat,Y_mat,X_prev,Y_prev));
    Frob_diff = Frob_diff./norm(X_prev,'fro')*norm(Y_prev,'fro');
    
%     progress = [progress; Frob_diff]
    if Frob_diff<tol,
        break;
    end
    
end
t=toc;
disp(['KGL timing = ' num2str(t) ' secs']);
end
