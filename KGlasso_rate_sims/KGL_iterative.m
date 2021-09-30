function [Frob_error_inv, Frob_error_cov, X_mat, Y_mat] = KGL_iterative(SCM,p,f,n,cx,cy,A0,B0,X0,Y0,N_iter,tol)
%KGL Function implements the iterative KGL algorithm and outputs the
%Frobenius norm errors.
%
% Last updated: Nov. 5, 2011
%

X_mat = eye(p); % starting point
Y_mat = eye(f);

% sparsify each factor
M = max(max(p,f),n);
lambda_initial = cy*sqrt(log(M)/(n*p));
lambda_final   = lambda_initial + cx*sqrt(log(M)/(n*f));
% lambda_x = cx*(1/sqrt(p)+1/sqrt(f))*sqrt(log(M)/n);
% lambda_y = cy*(1/sqrt(p)+1/sqrt(f))*sqrt(log(M)/n);

% progress = [];
for n1=1:N_iter,
    X_prev = X_mat;
    Y_prev = Y_mat;
    

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
        [B_hat_inv B_hat] = quicGlasso(T_X, lambda, lambda*10, 1e-3, 300);
%         B_hat_inv = glasso2(T_X, lambda);
%         B_hat = inv(B_hat_inv);
    else
        B_hat = T_X;
        B_hat_inv = inv(T_X);
    end
    inv_Y_mat = B_hat;
    Y_mat = B_hat_inv;
    
    
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
        [A_hat_inv A_hat] = quicGlasso(T_Y, lambda, lambda*10, 1e-3, 300);
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

    % check exit condition
    Frob_diff = sqrt(computeFrob(X_mat,Y_mat,X_prev,Y_prev));
    Frob_diff = Frob_diff/(norm(X_prev,'fro')*norm(Y_prev,'fro'));
    
%     progress = [progress; Frob_diff]
    if Frob_diff<tol,
        break;
    end
    
end

% figure, plot(progress), xlabel('Iteration');
% pause

% error computation
temp1 = (norm(X0,'fro')*norm(Y0,'fro'))^2;
temp2 = (norm(A0,'fro')*norm(B0,'fro'))^2;
Frob_error_inv = computeFrob(X_mat,Y_mat,X0,Y0)/temp1;
Frob_error_cov = computeFrob(inv_X_mat,inv_Y_mat,A0,B0)/temp2;
%     Frob_error_inv = norm(kron(X_mat,Y_mat)-kron(X0,Y0),'fro')^2/temp1
%     Frob_error_cov = norm(kron(inv_X_mat,inv_Y_mat)-kron(A0,B0),'fro')^2/temp2
%
%     pause


end

