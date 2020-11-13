function [Frob_error_ML_inv_final, Frob_error_ML_cov_final, A_mat, B_mat] = FF(SCM,p,f,n,niter, tol,A0,B0,X0,Y0)
%FF Function implements the flip-flop algorithm.
% 
% Last updated: Sept. 18, 2011
% 


% initial estimates of X, Y matrices
X_mat = eye(p);
Y_mat = eye(f);

for ni=1:niter,
    X_prev = X_mat;
    Y_prev = Y_mat;
        
    % compute temporary matrix T_X
    T_X = zeros(f,f);
    for i1=1:p,
        for i2=1:p,
        	T_X = T_X + X_mat(i2,i1)*SCM((i1-1)*f+1:i1*f,(i2-1)*f+1:i2*f);
        end
    end
    T_X = T_X/p;
    B_hat = T_X;
    B_hat_inv = inv(T_X);    
    Y_mat = B_hat_inv;
    
    % compute T_Y
    T_Y = zeros(p,p);
    for j1=1:f,
        for j2=1:f,
        	T_Y = T_Y + Y_mat(j2,j1)*SCM(j1:f:end,j2:f:end);
        end
    end
    T_Y = T_Y/f;
    A_hat = T_Y;
    A_hat_inv = inv(T_Y);
    X_mat = A_hat_inv;
    
    Frob_diff = sqrt(computeFrob(X_mat,Y_mat,X_prev,Y_prev));
    Frob_diff = Frob_diff/(norm(X_prev,'fro')*norm(Y_prev,'fro'));
    
    if Frob_diff<tol,
        break;
    end
    
end

    
temp1 = (norm(X0,'fro')*norm(Y0,'fro'))^2;
temp2 = (norm(A0,'fro')*norm(B0,'fro'))^2;
Frob_error_ML_inv_final = computeFrob(X_mat,Y_mat,X0,Y0)/temp1; % inverse
Frob_error_ML_cov_final = computeFrob(A_hat,B_hat,A0,B0)/temp2; % forward


A_mat = A_hat;
B_mat = B_hat;



end

