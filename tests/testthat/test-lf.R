library(romy)

test_that("check_lf", {
  N=1000
  M=30
  P_m=2
  P= P_m*M
  K=2
  beta=2.0
  X=matrix(rnorm(N*P),nrow=N,ncol=P)
  Y=matrix(rnorm(N*M),nrow=N,ncol=M)
  colnames(X)=paste0("X",1:ncol(X))
  colnames(Y)=paste0("Y",1:ncol(Y))

  inds_X=list()
  for(m in 1:M)
  {
		inds_X[[m]]=c(((m-1)*P_m+1):((m-1)*P_m+P_m))
  }
  Z1=matrix(0, nrow=N, ncol=K+1)
  F1=matrix(0, nrow=K+1, ncol=P)
			
  for(m in 1:M)
  {
			  alpha_m=rnorm(P_m, mean=0, sd=1)
			  F1[1,((m-1)*P_m+1):((m-1)*P_m+P_m)]=alpha_m
			  for(k in 1:K)
			  {
				beta_km=rnorm(P_m, mean=0, sd=1)
				F1[k+1,((m-1)*P_m+1):((m-1)*P_m+P_m)]=beta_km
			  }
  }
			
  Z1[,1]=1
  Z=matrix(rnorm(N*K, mean = 0, sd=1), nrow=N, ncol=K)
  Z1[,2:(K+1)]=Z
  F=F1[2:(K+1),]
			
			
  Gamma=Z1 %*% F1
			
			
			
  for(m in 1:M)
  {
			  Y[,m]=Y[,m]+beta*rowSums(as.matrix(X[,inds_X[[m]]]) * as.matrix(Gamma[,inds_X[[m]]]))
  }
  lf=romy_lf(Y=Y, X=X, inds_X=inds_X)
  
  expect_false(any(is.na(lf)))
  

})
