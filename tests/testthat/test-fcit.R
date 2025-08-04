library(romy)

test_that("check_fcit", {
  X=matrix(rnorm(1000*1),nrow=1000,ncol=1)
  Z=matrix(rnorm(1000*2),nrow=1000,ncol=2)
  Y=matrix(rnorm(1000*10),nrow=1000,ncol=10)
  colnames(X)=paste0("X",1:ncol(X))
  colnames(Y)=paste0("Y",1:ncol(Y))

  fcit=romy_fcit(Y=Y, X=X, Z=Z)
  
  expect_false(any(is.na(fcit)))
  

})