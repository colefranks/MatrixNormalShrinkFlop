function B =  GeminiB(X, rowpen)
  [p,f,n]=size(X);
  Gamma = zeros(p,p);
  for i=1:n
    Gamma = Gamma + X(:,:,i)*X(:,:,i)';
  end
  %this WILL penalize the diagonal
  
  glassoB = quicGlasso(Gamma,rowpen,rowpen*2,1e-3,300);
  
end

