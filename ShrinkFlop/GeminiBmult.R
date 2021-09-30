#' Estimate Row-Row Covariance Structure Using Gemini
#'
#' GeminiB estimates the row-row covariance, inverse covariance,
#' correlation, and inverse correlation matrices using Gemini.
#' For identifiability, the covariance factors A and B are scaled so
#' that A has trace m, where m is the number of columns of X,
#' A is the column-column covariance matrix, and B is the row-row
#' covariance matrix.
#'
#' @param kraus tuple of Data matrices, of dimensions n by m.
#' @param rowpen Glasso penalty parameter.
#' @param penalize.diagonal Logical value indicating whether to penalize the
#' off-diagonal entries of the correlation matrix.  Default is FALSE.
#' @return
#' \item{B.hat.inv}{estimated inverse covariance matrix.}
#' \item{A.hat.inv}{estimated inverse covariance matrix. }
#' @examples
#' n1 <- 5
#' n2 <- 5
#' n <- n1 + n2
#' m <- 20
#' X <- matrix(rnorm(n * m), nrow=n, ncol=m)
#' rowpen <- sqrt(log(m) / n)
#' out <- GeminiB(X, rowpen, penalize.diagonal=FALSE)
#' # Display the estimated correlation matrix rounded to two
#' # decimal places.
#' print(round(out$corr.B.hat, 2))
#' @export
GeminiBmult <- function (kraus, rowpen, penalize.diagonal=FALSE) {
  
  n = nrow(kraus[[1]])
  #print(n)
  m = ncol(kraus[[1]])
  #print(m)
  
  marginal1 = cpmap(kraus,diag(m))
  Gamma.hat.B <- stats::cov2cor(marginal1)
  #tic("entering glasso")
  glasso.B <- glasso::glasso(Gamma.hat.B, rowpen,
                             penalize.diagonal=penalize.diagonal)
  #toc()

  sd.row <- sqrt(diag(marginal1))
  inv.sd.row <- 1 / sd.row

  B.hat <- t(t(glasso.B$w * sd.row) * sd.row) / m
  B.hat.inv <- m * t(t(glasso.B$wi * inv.sd.row) * inv.sd.row)

  output <- list(corr.B.hat=glasso.B$w,
                 corr.B.hat.inv=glasso.B$wi,
                 B.hat=B.hat,
                 B.hat.inv=B.hat.inv)
}

GeminiBflop <- function (kraus, rowpenA, rowpenB, penalize.diagonal=FALSE) {
  niter = length(rowpenA)
  n = nrow(kraus[[1]])
  #print(n)
  m = ncol(kraus[[1]])
  
  r = length(kraus)
  #print(r)
  
  marginal1 = cpmap(kraus,diag(m))
  #marginal2 = cpdual(kraus,diag(n))
  
  sd.row <- sqrt(diag(marginal1))
  #sd.col <- sqrt(diag(marginal2))
  
  inv.sd.row <- 1 / sd.row
  #inv.sd.col <- 1 / sd.col
  
  Binv2 = diag(n)
  Ainv2 = diag(m)
  for (i in 1:(niter-1)){
    tic("marginal1")
    marginal1 = cpmap(kraus,solve(Ainv2))/r*m
    sd.row <- sqrt(diag(marginal1))
    inv.sd.row <- 1 / sd.row
    
    Gamma.hat.B <- diag(inv.sd.row)%*%marginal1%*%diag(inv.sd.row)
    glasso.B <- glasso::glasso(Gamma.hat.B, rowpenA[i],
                               penalize.diagonal=penalize.diagonal)
    
    #B.hat <- t(t(glasso.B$w * sd.row) * sd.row) / m
    Binv2 <- m * t(t(glasso.B$wi * inv.sd.row) * inv.sd.row)
    
    toc()
    
    tic("marginal2")
    marginal2 = cpdual(kraus,solve(Binv2))/r*n
    sd.col <- sqrt(diag(marginal2))
    inv.sd.col <- 1 / sd.col
    
    Gamma.hat.A <- diag(inv.sd.col)%*%marginal2%*%diag(inv.sd.col)
    glasso.A <- glasso::glasso(Gamma.hat.A, rowpenB[i],
                               penalize.diagonal=penalize.diagonal)
    Ainv2 <- n * t(t(glasso.A$wi * inv.sd.col) * inv.sd.col)
    toc()
  }
  marginal1 = cpmap(kraus,solve(Ainv2))/r*m
  sd.row <- sqrt(diag(marginal1))
  inv.sd.row <- 1 / sd.row
  
  Gamma.hat.B <- diag(inv.sd.row)%*%marginal1%*%diag(inv.sd.row)
  glasso.B <- glasso::glasso(Gamma.hat.B, rowpenA[i],
                             penalize.diagonal=penalize.diagonal)
  
  #B.hat <- t(t(glasso.B$w * sd.row) * sd.row) / m
  Binv2 <- m * t(t(glasso.B$wi * inv.sd.row) * inv.sd.row)
  
  output <- list(A.hat.inv2=Ainv2,
                 B.hat.inv2=Binv2)
}
