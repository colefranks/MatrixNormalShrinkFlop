n = 2
m = 3
A = matrix(rnorm(n * m), nrow=n, ncol=m)
B = matrix(rnorm(n * m), nrow=n, ncol=m)
krauses = list(A,B)

sum(krauses %*% t(krauses))
krauses %*% lapply(krauses, t)

cpmap(t(A))
krauses
lapply(krauses, t)

y=0
for(x in krauses){
  y<-y + x%*% t(x)
}
y  
dim(krauses[1])
dim(A)[2]

sinkstep(krauses, diag(m))
solve(A %*% t(A))/2
C = krauses[[1]]
C
A
dim(krauses[[1]])[[1]]
cpdual(krauses,diag(3))
krauses
lapply(krauses,t)
A
t(A)%*%diag(2)
tr(A%*%t(A))
A%*%t(A)
sinkhorn(krauses,verbose=TRUE)
#question - is frob also a rep? tr g^Tg is a norm, and to any power I suppose it 
#would be too, though sinkhorn would be less clear. g^Tg g^Tg on the other hand, idk.
# is an inner product - <X, X> under conjugation action, which I guess is also a rep.
regsinkhorn(krauses,verbose=TRUE,reg=0.01)

n = 10
m = 20

A = .2*diag(n)
for (term in c(1:s)){
  #print(term)
  x=sample(1:n,2,replace=F)
  A[x[1],x[1]]=A[x[1],x[1]]+1
  A[x[1],x[2]]=A[x[1],x[2]]-1
  A[x[2],x[1]]=A[x[2],x[1]]-1
  A[x[2],x[2]]=A[x[2],x[2]]+1
}
A1 = sqrtm(solve(A))

B = .2*diag(m)
for (term in c(1:s)){
  #print(term)
  x=sample(1:n,2,replace=F)
  B[x[1],x[1]]=B[x[1],x[1]]+1
  B[x[1],x[2]]=B[x[1],x[2]]-1
  B[x[2],x[1]]=B[x[2],x[1]]-1
  B[x[2],x[2]]=B[x[2],x[2]]+1
}
B1 = sqrtm(solve(B))
X1 <- matrix(rnorm(n * m), nrow=n, ncol=m)
X2 <-matrix(rnorm(n * m), nrow=n, ncol=m)
X1<- A1 %*% X1 %*% B1
X2<- A1 %*% X2 %*% B1

#rowpen.list <- sqrt(log(m) / n) * c(1, 0.5, 0.1)
#out <- GeminiBPath(X, rowpen.list, penalize.diagonal=FALSE)

rowpen <- sqrt(log(m) / n) 
rowpen2 <-rowpen*c(.1,1)
out <- GeminiBmult(list(X1), rowpen, penalize.diagonal=FALSE)
out2 <- GeminiBflop(list(X1), rowpen2, penalize.diagonal=FALSE)
out3<-GeminiB(X1,rowpen,penalize.diagonal=FALSE)
E = n*A/tr(A)
C = out$B.hat.inv
C = n*C/tr(C)
D = regsinkhorn(list(X1,X2),reg=regu,maxit=3)[[1]]
D = n*D/tr(D)
G = out2$B.hat.inv
G = n*G/tr(G)
H = out3$B.hat.inv
H = n*H/tr(H)

round(E,2)
round(C,2)
round(D,2)
round(G,2)
round(H,2)

thresh(.1*diag(3),.2)
GeminiB(c(X1,X2), rowpen, penalize.diagonal=FALSE)
out