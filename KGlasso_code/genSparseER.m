function [A0,X0] = genSparseER(p)
%GENSPARSEER Function generates sparse matrix using Erdos Renyi model.

probA = 0.8;

Ainv = rand(p,p);
Ainv = Ainv>probA;
Ainv = (Ainv+Ainv.')/2;
Ainv(Ainv==0.5)=0;
% Ainv(linspace(1,numel(Ainv),length(Ainv))) = 1; % replace diagonal with 1's
% Binv(linspace(1,numel(Binv),length(Binv))) = 1; % replace diagonal with 1's
% make positive definite
% epsilon = 1;
epsilon = 0.5;
rho = epsilon - min(eig(Ainv));
Ainv = Ainv + rho*eye(p);
X0 = Ainv./norm(Ainv,2);
A0 = inv(X0);

end

