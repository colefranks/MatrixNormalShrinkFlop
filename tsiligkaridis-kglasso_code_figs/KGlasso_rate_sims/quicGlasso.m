function [X,W,numIter] = quicGlasso(SCM, lambda, tol, dgap_tol, maxIter)
%QUICGLASSO MATLAB implementation of QUIC (Hsieh,et al. NIPS 2011).

p = size(SCM,1);
Xprev = eye(p); % initial guess for inverse covariance
Wprev = eye(p); % forward covariance
beta = 0.5;
sigma = 0.25;
epsilon = 1e-2; % inner stopping tolerance
for m=1:maxIter,
    D = zeros(p,p);
    U = zeros(p,p); % slack variable
    grad_g = -Wprev + SCM;
%     figure(21), imagesc(grad_g), axis square, pause;
    Sfree_I = [];
    Sfree_J = [];
    for i1=1:p,
        for i2=1:p,
            if (abs(grad_g(i1,i2))>=lambda-epsilon) || (Xprev(i1,i2)~=0),
                Sfree_I = [Sfree_I; i1];
                Sfree_J = [Sfree_J; i2];
            end
        end
    end
%     size(Sfree_J)
%     figure, hist(SfreeI)
%     figure, hist(Sfree_I), pause
%     [Sfree_I Sfree_J]', pause
    
    diffD = 0;
    normD = 0;
    Dnew = zeros(p,p);
    Unew = zeros(p,p);
    
    while(1),
        Dprev = D;
        Uprev = U;
        
        for ii=1:length(Sfree_I),
            i = Sfree_I(ii);
            j = Sfree_J(ii);
            
            if i<=j, % ensure symmetry of D
                a = Wprev(i,j)^2 + Wprev(i,i)*Wprev(j,j);
    %             if i~=j,
    %                 a = a + Wprev(i,i)*Wprev(j,j);
    %             end
                b = SCM(i,j) - Wprev(i,j) + Wprev(:,i)'*U(:,j);
    %             b = SCM(i,j) - Wprev(i,j) +  Wprev(:,i)'*D*Wprev(:,j);

                c = Xprev(i,j) + D(i,j);
%                 mu = -c + softThres(c-b/a,lambda/a);
                mu = -c + sign(c-b/a)*max(abs(c-b/a)-lambda/a,0);
                D(i,j) = D(i,j) + mu; % update D symmetrically
                D(j,i) = D(j,i) + mu;


%                 normD = normD - abs(D(i,j));

%                 l = lambda/a;
%                 f = b/a;
%                 if c>f,
%                     mu = -f-l;
%                     if c+mu<0,
%                         mu = -c;
%                         D(i,j) = -Xprev(i,j); % direct
%                     else
%                         D(i,j) = D(i,j) + mu;
%                     end
%                 else
%                     mu = -f+l;
%                     if c+mu>0,
%                         mu = -c;
%                         D(i,j) = -Xprev(i,j);
%                         
%                     else
%                         D(i,j) = D(i,j) + mu;
%                     end
%                 end
%                 D(j,i) = D(i,j);

                diffD = diffD + abs(mu);
                normD = normD + abs(D(i,j));

                U(i,:) = U(i,:) + mu*Wprev(j,:);
                U(j,:) = U(j,:) + mu*Wprev(i,:);
            end
        end
                
%         D = D + D';
%          figure(20), imagesc(D), axis square, colormap gray, colormap(1-colormap), title('Matrix D'), pause;
   
%         errD = norm(Dprev-D,'fro');
%         errU = norm(Uprev-U,'fro');
%         [errD errU], pause
%         if errD<tol,
%             break;
%         end

%         [diffD normD*tol], pause
        if diffD < normD*tol,
            break;
        end
        

%         disp('finished iter'), pause;
    end
    
    
%     % D should be symmetric - but is not!!
%     disp(['Newton iteration = ' num2str(m)]);
%     figure(20), imagesc(D), axis square, colormap gray, colormap(1-colormap), title('Matrix D'), pause;
%     [min(abs(D(:))) max(abs(D(:)))]

    f_prev = -log(det(Xprev)) + trace(SCM*Xprev) + lambda*norm(Xprev(:),1);
    power = -1;
    while(1),
        power = power + 1;
        alpha = beta^power;
%         min(eig(Xprev))
%         min(eig(Xprev+alpha*D))
%         pause
%         D
%         power, pause
        
        if min(eig(Xprev+alpha*D))<=0,
            continue;
        end
%         eig(Xprev+alpha*D), pause
        L = chol(Xprev+alpha*D,'lower');
        f_current = -2*log(det(L)) + trace(SCM*L*L') + lambda*norm(Xprev+alpha*D,1);
        if f_current<=f_prev + alpha*sigma*(trace(grad_g*D)+lambda*norm(Xprev(:)+alpha*D(:),1)-lambda*norm(Xprev(:),1)),
            break;
        end        
    end
    
    % update estimates
    X = Xprev + alpha*D;
    W = inv(L')*inv(L);
%     W = inv(Xprev+alpha*D);

    dgap = trace(SCM*X)+lambda*norm(X(:),1)-p;
    if abs(dgap)<dgap_tol, % exit based on duality gap
        break;
    end
    
    % prepare for next iteration
    Xprev = X;
    Wprev = W;
end

% number of Newton iterations until exit
numIter = m;


end % close function

