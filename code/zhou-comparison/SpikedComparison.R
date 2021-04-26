library(jointMeanCov)
library(Matrix)
library(expm)
library(matrixcalc)
library(ggplot2)
library(tictoc)


#normalize determinant

DetNorm <-function(A){
  n = dim(A)[1]
  return(A*det(A)^(-1/n))
}

#define geodesic distance
GeodesicDistance <-function(P,Q){
  root = solve(sqrtm(Q))
  P = root %*% P %*% root
  return(norm(logm(P),"F"))
}

#test geodesic distance
A = diag(c(1,2,3))
A = expm(A)
GeodesicDistance(diag(c(1,1,1)), A)

SpikedComparison <-function(SmallDimension,BigDimension,NumSamples,NumInstances=2,TrialsPerInstance=2,spike = 10, RegMin = -5, RegMax =5, RegStride=1){
  
  Regularizers = exp(RegStride*c(floor(RegMin/RegStride):floor(RegMax/RegStride)))
  print(Regularizers)
  results = data.frame(Regularizer = Regularizers,Gemini=0,RegSink = 0,Trivial=0)
  GeoResults = data.frame(Regularizer = Regularizers,Gemini=0,RegSink = 0,Trivial=0)
  OpResults = data.frame(Regularizer = Regularizers,Gemini=0,RegSink = 0,Trivial=0)
  
  for (j in 1:length(Regularizers)){
    
    #for each value of the regularizer we will compute the error of the estimator on several ground truth covariances

    for (l in 1:NumInstances){ 
  
      #define the covariances
      A = diag(SmallDimension)
      v = matrix(rnorm(SmallDimension),nrow = SmallDimension, ncol = 1)
      A = A + spike*v%*%t(v)
      rootA = sqrtm(A)
      invA = solve(A)
      invA = invA/tr(invA)
      normA = norm(invA,"F")^2
    
      B = diag(BigDimension)
      v = matrix(rnorm(BigDimension),nrow = BigDimension, ncol = 1)
      B = B + spike*v%*%t(v)
      rootB = sqrtm(B)
      invB = solve(B)
    
      
      #for each covariance we will compute the estimate with fresh samples from the model a few times
      for(i in 1:TrialsPerInstance){
        
        #create the fresh samples
        X = list(0*c(1:NumSamples))
        
        for(i in 1:NumSamples){
          
          #create the data
          X[[i]] = rootA%*%matrix(rnorm(SmallDimension*BigDimension,sd=1), nrow=SmallDimension, ncol=BigDimension)%*%rootB
          
          
        }
        
        #trivial frobenius error
        TrivialError = norm(invA - diag(SmallDimension)/SmallDimension,"F")^2 / normA
        results$Trivial[j] = results$Trivial[j]+ TrivialError/(NumInstances*TrialsPerInstance)
        
        #operator error
        OpTrivialError = norm(invA - diag(SmallDimension)/SmallDimension,"I")^2
        OpResults$Trivial[j] = OpResults$Trivial[j]+ OpTrivialError/(NumInstances*TrialsPerInstance)
        
        #geodesic trivial
        GeoTrivialError = GeodesicDistance(DetNorm(invA), diag(SmallDimension))^2
        GeoResults$Trivial[j] = GeoResults$Trivial[j] + GeoTrivialError/(NumInstances*TrialsPerInstance)
        
        
        #gemini
        #set gemini penalties
        GeminiRegFactor = 1*sqrt(log(SmallDimension)/BigDimension)
        tic("geminiB")
        out <- GeminiBmult(X, GeminiRegFactor*Regularizers[j], penalize.diagonal=FALSE)
        toc()
        GeminiEstimate = out$B.hat.inv
        #calculate frobenius error
        GeminiError = norm(invA - (GeminiEstimate/tr(GeminiEstimate)),"F")^2 / normA
        results$Gemini[j] = results$Gemini[j]+ GeminiError/(NumInstances*TrialsPerInstance)
        
        #operator error
        OpGeminiError = norm(invA - (GeminiEstimate/tr(GeminiEstimate)),"I")^2
        OpResults$Gemini[j] = OpResults$Gemini[j]+ OpGeminiError/(NumInstances*TrialsPerInstance)
        
        #calculate geodesic error
        GeoGeminiError = GeodesicDistance(DetNorm(invA), DetNorm(GeminiEstimate))^2
        GeoResults$Gemini[j] = GeoResults$Gemini[j] + GeoGeminiError/(NumInstances*TrialsPerInstance)
        
        #regularized Sinkhorn
        #set sinkhorn penalties
        RegSinkFactor = 10
        tic("regularized sinkhorn")
        RegSinkEstimate = regsinkhorn(X,tol = .1*Regularizers[j], reg = Regularizers[j]*RegSinkFactor)[[1]]
        toc()
        
        #calculate frobenius
        RegSinkError = norm(invA - (RegSinkEstimate/tr(RegSinkEstimate)),"F")^2 / normA
        results$RegSink[j] = results$RegSink[j]+ RegSinkError/(NumInstances*TrialsPerInstance)
        
        #operator error
        OpSinkError = norm(invA - (RegSinkEstimate/tr(RegSinkEstimate)),"I")^2
        OpResults$RegSink[j] = OpResults$RegSink[j]+ OpSinkError/(NumInstances*TrialsPerInstance)
        
        #geodesic distance
        GeoSinkError = GeodesicDistance(DetNorm(invA), DetNorm(RegSinkEstimate))^2
        GeoResults$RegSink[j] = GeoResults$RegSink[j] + GeoSinkError/(NumInstances*TrialsPerInstance)
        
        

      }
    }
  }
  output<- list(Results = results,OpResults = OpResults, GeoResults = GeoResults)
}

ans = SpikedComparison(25,50,1,RegMin=-10,NumInstances =2, TrialsPerInstance  = 2, RegStride = 1)
SimplePlot(ans$GeoResults)

