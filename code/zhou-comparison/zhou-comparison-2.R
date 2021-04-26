library(jointMeanCov)
library(Matrix)
library(expm)
library(matrixcalc)
library(ggplot2)
#library(devtools)
library(tictoc)
#have had some trouble with this hanging. Seems to happen at random.

#n1 <- 5
#n2 <- 5
n <- 10
m <- 20
s=4
k = 3
ntrials = 2
ninstances = 2
regus = exp(1*c(-10:0))
errors = matrix(0,nrow = 4,ncol = length(regus))
results = data.frame("Regularizer"=regus)
for (l in 1:ninstances){

  A = .05*diag(n)
  for (term in 1:s*n/10){
    #print(term)
    x=sample(1:n,2,replace=F)
    A[x[1],x[1]]=A[x[1],x[1]]+1
    A[x[1],x[2]]=A[x[1],x[2]]-1
    A[x[2],x[1]]=A[x[2],x[1]]-1
    A[x[2],x[2]]=A[x[2],x[2]]+1
  }
#add big dense piece spike
  v = matrix(0,n)
  v[1:n/2] = 1

  A = A + 50*v%*%t(v)

#A[1:20,1:20]


  A1 = sqrtm(solve(A))

  B = .05*diag(m)
  for (term in 1:s*m/10){
    #print(term)
    x=sample(1:m,2,replace=F)
    B[x[1],x[1]]=B[x[1],x[1]]+1
    B[x[1],x[2]]=B[x[1],x[2]]-1
    B[x[2],x[1]]=B[x[2],x[1]]-1
    B[x[2],x[2]]=B[x[2],x[2]]+1
  }

#v = matrix(0,50)
#v[1:45] = 1

#add dense thing to other side
  v = matrix(rnorm(m*(m-1)),nrow = m, ncol = m-1)
  B = B + 1000*v%*%t(v)


#make other side ill-conditioned
#v = matrix(0,m)
#v[1:(m-5)]=1
#B = B + 200*diag(c(v))

  B1 = sqrtm(solve(B))


  for (j in 1:length(regus)){
    regu = regus[j]
    print(c("REG:",regu,"RUN",j))
    rmse1 = 0
    rmse2 = 0
    rmse3 = 0
    rmse4 = 0
    for(i in 1:ntrials){
      print(c("trial",i))
  
 
      X = list(0*c(1:k))
      for(i in 1:k){
        X[[i]] = A1%*%matrix(rnorm(n*m,sd=1), nrow=n, ncol=m)%*%B1
        #print(length(X))
      }
      print(c(length(X),"X length"))
      # X<- A1 %*% X %*% B1
    
      #rowpen.list <- sqrt(log(m) / n) * c(1, 0.5, 0.1)
      #out <- GeminiBPath(X, rowpen.list, penalize.diagonal=FALSE)
      #tic("sinkhorn")
      max = 10*(1+ log(regu/min(regus)))
      #print(max)
    
    

    
      rowpen <- sqrt(log(m) / n) 
      rowpenA<-rowpen*c(200,200)
      rowpenB<-rowpen*c(2000000,2000000)
      tic("geminiB")
      out <- GeminiBmult(X, rowpen*regu, penalize.diagonal=FALSE)
      toc()
    
      tic("sinkhorn")
      #D = threshsinkhorn(X,thr=0,tol=.1*regu,reg=1)[[1]]
      #D = thresh(n*D/tr(D),.00)
      #D = regsinkhorn(X,tol=.1*regu,reg=regu/500)[[1]]
      #D = regsinkhorn(X,tol = .1*regu,reg = 10*regu, verbose=TRUE)[[1]]
      #D = n*D/tr(D)
      #D = thresh(n*D/tr(D),.00)
      #D = sinkhornlast(X,tol = .1*regu,thr = 0.1, reg = 10*regu, verbose=TRUE)[[1]]
      #D = glassosinkhorn(X,tol = .1*regu, reg = 10*regu, glassoreg = rowpen*regu/5, verbose=TRUE)[[1]]
      #D = glassosinkhorn(X,tol = .1*regu, reg = regu/1000, glassoreg = rowpen*regu/5, verbose=TRUE)[[1]]
      D = regsinkhorn(X,tol = .1*regu, reg = regu/4000)[[1]]
      D = n*D/tr(D)
      toc()

      tic("geminiBflop")
      #out2 <-GeminiBflop(X,rowpenA*regu/5, rowpenB*regu/5,penalize.diagonal=FALSE)
      toc()
    
      E = n*A/tr(A)
      C = out$B.hat.inv
      C = n*C/tr(C)
    
      #H = out2$B.hat.inv2
      H = matrix(0,nrow=n,ncol=n)
      #H = n*H/tr(H)
    
      if(j==1){
        E1 = E
        C1 = C
        D1 = D
        H1= H
      }
      #H = out2$B.hat.inv2
      #H = n*H/tr(H)
  
      #K = regsinkhorn(list(X),reg=1,maxit=2)[[1]]
      #K = thresh(n*K/tr(K),0)
    
    
    
      norm = norm(E,"F")^2
      rmse1= rmse1 + norm(C- E, "F")^2/norm
      rmse2 = rmse2 + norm(D-E, "F")^2/norm
      rmse3 = rmse3 + norm(diag(n)-E, "F")^2/norm
      rmse4 = rmse4 + norm(H-E, "F")^2/norm
    
    }

    errors[,j]=errors[,j] + c(rmse1/(ntrials*ninstances),
               rmse2/(ntrials*ninstances),
               rmse3/(ntrials*ninstances),
               rmse4/(ntrials*ninstances)
               )
  #print(errors)
  #print(errors)
  #rmse1/ntrials #theirs
  #round(E,2)[1:5,1:5]
  #rmse2/ntrials #ours
  #round(C,2)[1:5,1:5]
  #rmse3/ntrials #trivial thing - how far it is from identity
  #round(D,2)[1:5,1:5]
}
}
results = cbind(results,"Regularized Sinkhorn"=errors[2,])
results = cbind(results,"Zhou"=errors[1,])

SimplePlot(results)

errors
#it does seem that frobenius is marginally better. 


#C
#round(out$B.hat.inv,2)
# ok this thing is pretty impressive at nailing the graph structure, but not the weights. 
#the frobenius thing seems really bad.


#next need to look at zhou's flip flop. I think it should just be alternate applications of 
#geminiB. 

#should also, I suppose, compare with allen-tang. 
#plot(c(1:10))
#par(new=FALSE)
#plot(2*c(1:10), add=TRUE)
#points(c(1:10))
#plot(sin,-pi, 4*pi, col = "red")
#plot(cos,-pi, 4*pi, col = "blue", add = TRUE)
#points(c(1:10))
#plot(1:10, 0:3, type = "n")
#plot(y = 0:length(regus),x = 0:length(regus),type="n")
#plot(1.0:10.0,0.0:2.0)
ran = max(errors) - min(errors)
plot(regus,log(min(errors) + (ran/length(regus))*c(0:(length(regus)-1))),type="n", log = "x")
#plot(regus, errors[1,])
#points(10*errors[1,],regus)
points(regus, log(errors[1,]))
points(regus, log(errors[2,]), col='red')
points(regus, log(errors[3,]), col='blue')
#points(regus, log(errors[4,]), col='green')
round(E1,2)[1:20,1:20]
round(D1,2)[1:20,1:20]
round(C1,2)[1:20,1:20]

#yay! in certain situations, it's better to just do glasso at the end. 