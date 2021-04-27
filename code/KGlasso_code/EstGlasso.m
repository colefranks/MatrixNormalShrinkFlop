function [Frob_err_inv, Frob_err_cov] = EstGlasso(SCM,lambdaGlasso,A0,B0,X0,Y0)
%ESTGLASSO Function for naive Glasso.

[Cov_est_inv Cov_est] = quicGlasso(SCM,lambdaGlasso,lambdaGlasso*2,1e-3,300);
Frob_err_inv = norm(Cov_est_inv-kron(X0,Y0),'fro')^2/norm(kron(X0,Y0),'fro')^2;
Frob_err_cov = norm(Cov_est-kron(A0,B0),'fro')^2/norm(kron(A0,B0),'fro')^2;

end

