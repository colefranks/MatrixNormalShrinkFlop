#
#completely positive map
cpmap<-function(kraus, current){
  y=0
  for(x in kraus){
    #print(dim(x))
    y<-y + x%*% current %*% t(x)
  }
  y
}

#dual of completely positive map
cpdual<-function(kraus, current){
  cpmap(lapply(kraus, t), current)
}

#trace -_-
tr<-function(mat){
  x = 0
  n = dim(mat)[1]
  for(i in 1:n){
    x = x+mat[i,i]
  }
  x
}

#one sinkhorn step
sinkstep<-function(kraus, current){
  #get dimensions
  m = dim(kraus[[1]])[1]
  #print(m)
  n = dim(kraus[[1]])[2]
  #print(n)
  #compute left scaling
  left = solve(cpmap(kraus,current))/m
  #print(left)
  #compute right scaling
  #print(cpdual(kraus,left))
  right = solve(cpdual(kraus, left))/n
  #output the pair of them
  list(left,right)
}

#sinkhorn's algorithm
sinkhorn<-function(kraus, tol=10**(-6),maxit=100,verbose=FALSE){
  m = dim(kraus[[1]])[1]
  n = dim(kraus[[1]])[2]
  #initialize number of iterations to zero
  it = 0
  #start with the identity
  left = diag(m)
  right = diag(n)
  #initialize eps to 1, something large
  eps = 1
  #while not scaled
  while(eps > tol & it < maxit){
    #measure progress
    eps = tr(matrix.power((left %*% cpmap(kraus,right)) - diag(m)/m, 2))
    if(verbose){
      print(eps)
    }
    #do a sinkhorn step
    sink = sinkstep(kraus,right)
    left = sink[[1]]
    right = sink[[2]]

    #increase iterations
    it =it+1
  }
  list(left,right)
}

#idea for regsink: add some of the map X \mapsto \tr(X) I.
#easiest way seems to be to define regcp and regcpdual.

#regularized completely positive map
regcpmap<-function(kraus, current, reg=0.1){
  y=0
  
  n = dim(kraus[[1]])[1]
  for(x in kraus){
    y<-y + x%*% current %*% t(x)
  }
  y + reg*tr(current)*diag(n)
}

#dual of regularized completely positive map
regcpdual<-function(kraus, current, reg=0.1){
  regcpmap(lapply(kraus, t), current, reg)
}

#one regularized sinkhorn step
regsinkstep<-function(kraus, current, reg){
  #get dimensions
  m = dim(kraus[[1]])[1]
  #print(m)
  n = dim(kraus[[1]])[2]
  
  numsamp = length(kraus)
  tracerho = tr(cpmap(kraus,diag(n)))
  
    
  #print(n)
  #compute left scaling
  #print(round(regcpmap(kraus,current,reg),2))
  left = solve((1 - reg)*cpmap(kraus,current)/(numsamp*m*n) + reg*tracerho*tr(current)*diag(m)/(m*n))/m
  #print(round(m*left,2))
  #print(left)
  #compute right scaling
  #print(cpdual(kraus,left))
  right = solve((1 - reg)*cpdual(kraus,left)/(numsamp*m*n) + reg*tracerho*tr(left)*diag(n)/(m*n))/n
  #output the pair of them
  list(left,right)
}

#do for tensors also?



threshsinkstep<-function(kraus, current, thr, reg){
  #print(thr)
  #get dimensions
  m = dim(kraus[[1]])[1]
  #print(m)
  n = dim(kraus[[1]])[2]
  #print(n)
  #compute left scaling
  #tic("cpmap")
  marginal=regcpmap(kraus,current,reg)
  #toc()
  
  #print(round(marginal,2))
  #print(round(solve(marginal),2))
  #print("marginal before funnies")
  #print(solve(marginal)[1:3,1:3])
  marginal = solve(marginal)
  #print(diag(marginal))
  sds = diag(diag(marginal)^(-1/2))
  #print("left")
  #print(sds)
  marginal = sds%*%marginal%*%sds
  #print(round(marginal,2))
  #print(marginal)
  #print("after 2 corr")
  #tic("thresholding")
  left = solve(sds)%*%psdthresh(marginal,thr)%*%solve(sds)/m
  #toc()
  #print(left)
  #print("left")
  #print(round(m*left,2))
  #print("marginal after funnies")
  #print(left)
  #compute right scaling
  #print(cpdual(kraus,left))
  
  marginal=regcpdual(kraus,left,reg)
  
  #print(matrix.rank(marginal))
  #print(solve(marginal)[1:3,1:3])
  marginal = solve(marginal)
  #print(marginal)
  #print(diag(marginal))
  sds = diag(diag(marginal)^(-1/2))
  #print(sds)
  #print("right")
  marginal = sds%*%marginal%*%sds
  #print(marginal)
  right = (solve(sds))%*%psdthresh(marginal,thr)%*%(solve(sds))/n
  #output the pair of them
  list(left,right)
}


#regularized sinkhorn

regsinkhorn<-function(kraus, tol=10**(-6),maxit=100,verbose=FALSE,reg=0.1){
  m = dim(kraus[[1]])[1]
  n = dim(kraus[[1]])[2]
  numsamp = length(kraus)
  tracerho = tr(cpmap(kraus,diag(n)))
  #initialize number of iterations to zero
  it = 0
  #start with the identity
  left = diag(m)
  right = diag(n)
  #initialize eps to 1, something large
  eps = 1
  #while not scaled
  reg1 = (2/pi)*atan(reg)
  print(c("regsink reg",reg1))
  while(eps > tol & it < maxit){
    
    #measure progress
    eps = tr(matrix.power(left %*% ((1 - reg1)*cpmap(kraus,right)/(m*n*numsamp) + reg1*tr(right)*tracerho*diag(m)/(m*n)) - diag(m)/m, 2))
    
    #do a sinkhorn step
    left = regsinkstep(kraus,right,reg1)[[1]]
    right = regsinkstep(kraus,right,reg1)[[2]]
    
    if(verbose){
      print(eps)
    }
    #increase iterations
    it =it+1
  }
  list(left,right)
}

threshsinkhorn<-function(kraus, tol=10**(-2),maxit=20,verbose=FALSE,thr=0.5,reg=1){
  m = dim(kraus[[1]])[1]
  n = dim(kraus[[1]])[2]
  #initialize number of iterations to zero
  it = 0
  #start with the identity
  left = diag(m)
  right = diag(n)
  #initialize eps to 1, something large
  eps = 1
  #while not scaled
  while(eps > tol & it < maxit){
    #do a sinkhorn step
    eps = tr(matrix.power(left %*% regcpmap(kraus,right,reg) - diag(m)/m, 2))
    if(verbose){
      print(eps)
    }
    left = threshsinkstep(kraus,right,thr,reg)[[1]]
    #print(left[1:4,1:4])
    right = threshsinkstep(kraus,right,thr,reg)[[2]]
    #print(right[1:4,1:4])
    
    #measure progress
    
    #increase iterations
    it =it+1
  }
  list(left,right)
}

sinkhornlast<-function(kraus, tol=10**(-6),maxit=100,verbose=FALSE,reg=0.1,thr=0.1){
  m = dim(kraus[[1]])[1]
  n = dim(kraus[[1]])[2]
  #initialize number of iterations to zero
  it = 0
  #start with the identity
  left = diag(m)
  right = diag(n)
  #initialize eps to 1, something large
  eps = 1
  #while not scaled
  while(eps > tol & it < maxit){
    
    #measure progress
    eps = tr(matrix.power(left %*% regcpmap(kraus,right,reg) - diag(m)/m, 2))
    
    #do a sinkhorn step
    left = regsinkstep(kraus,right,reg)[[1]]
    right = regsinkstep(kraus,right,reg)[[2]]
    
    if(verbose){
      print(eps)
    }
    #increase iterations
    it =it+1
  }
  #only does left right now
  sds = diag(diag(left)^(-1/2))
  invsds = diag(diag(left))
  #print("left")
  #print(sds)
  left = sds%*%left%*%sds
  #print(round(marginal,2))
  #print(marginal)
  #print("after 2 corr")
  #tic("thresholding")
  left = invsds%*%psdthresh(left,thr)%*%invsds/m
  list(left,right)
}

#how's glasso-ing the last step do?


glassosinkhorn<-function(kraus, tol=10**(-6),maxit=100,verbose=FALSE,reg=0.1,glassoreg=1){
  m = dim(kraus[[1]])[1]
  n = dim(kraus[[1]])[2]
  #initialize number of iterations to zero
  it = 0
  #start with the identity
  left = diag(m)
  right = diag(n)
  #initialize eps to 1, something large
  eps = 1
  #while not scaled
  while(eps > tol & it < maxit){
    
    #measure progress
    eps = tr(matrix.power(left %*% regcpmap(kraus,right,reg) - diag(m)/m, 2))
    
    #do a sinkhorn step
    left = regsinkstep(kraus,right,reg)[[1]]
    right = regsinkstep(kraus,right,reg)[[2]]
    
    if(verbose){
      print(eps)
    }
    #increase iterations
    it =it+1
  }
  #only fixes left
  marginal1 = solve(left)
  sd.row <- sqrt(diag(marginal1))
  inv.sd.row <- 1 / sd.row
  
  Gamma.hat.B <- diag(inv.sd.row)%*%marginal1%*%diag(inv.sd.row)
  glasso.B <- glasso::glasso(Gamma.hat.B, glassoreg,
                             penalize.diagonal=FALSE)
  
  #B.hat <- t(t(glasso.B$w * sd.row) * sd.row) / m
  left <- m * t(t(glasso.B$wi * inv.sd.row) * inv.sd.row)
  list(left,right)
}


thresh<-function(x,tol){
  m = dim(x)[1]
  n = dim(x)[2]
  for(i in 1:n){
    for(j in 1:m){
      #print(x[i,j])
      if(abs(x[i,j])<tol){
        x[i,j]=0
      }
    }
  }
  return(x)
}

psdthresh<-function(x,tol){
  m = dim(x)[1]
  n = dim(x)[2]
  for(i in 2:n){
    for(j in 1:(i-1)){
      #print(x[i,j])
      if(abs(x[i,j])<tol){
        y = abs(x[i,j])
        x[i,j]=0
        x[i,i] = x[i,i]+y
        x[j,j] = x[j,j]+y
        x[j,i]=0
      }
    }
  }
  return(x)
}






zhouflop<-function(X,reg=what){
  
}



