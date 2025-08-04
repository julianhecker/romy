library(romy)

test_that("check_lf", {
  M=30
  P_m=2
  P= P_m*M
  X=matrix(rnorm(1000*P),nrow=1000,ncol=P)
  Y=matrix(rnorm(1000*M),nrow=1000,ncol=M)
  colnames(X)=paste0("X",1:ncol(X))
  colnames(Y)=paste0("Y",1:ncol(Y))

  inds_X=list()
  for(m in 1:M)
  {
		inds_X[[m]]=c(((m-1)*P_m+1):((m-1)*P_m+P_m))
  }
  lf=romy_lf(Y=Y, X=X, inds_X=inds_X)
  
  expect_false(any(is.na(lf)))
  

})
