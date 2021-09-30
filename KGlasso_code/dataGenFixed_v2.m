function [SCM] = dataGenFixed_v2(n,A0,B0)
% Function generates multivariate normal data based on given covariance
% matrix.

% global SCM;

p=size(A0,1);
f=size(B0,1);

% % whos
% % clear SCM;
% % Sigma_0 = kron(A0,B0); % stub
% 
% % generate data from multivariate normal
% mean0 = zeros(n,p*f); % 0-mean data
% data = mvnrnd(mean0,kron(A0,B0)).';

% alternative way to generate
data = [];
A0_sqrt = sqrtm(A0);
B0_sqrt = sqrtm(B0);
for m=1:n,
    W = randn(f,p);
    X = B0_sqrt*W*A0_sqrt;
    data = [data X(:)];
end
clear W X A0_sqrt B0_sqrt;
% clear Sigma_0;
% data(:,1:3), pause

% size(data_out) % pf x n
% SCM = cov(data.');
SCM = zeros(p*f,p*f);
for m=1:n,
    SCM = SCM + data(:,m)*data(:,m)';
end
SCM = SCM/n;


% SCM, pause

