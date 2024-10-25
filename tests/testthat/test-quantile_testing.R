library(romy)

test_that("check_quantile_testing", {
  X=matrix(rnorm(1000*3),nrow=1000,ncol=3)
  Z=matrix(rnorm(1000*3),nrow=1000,ncol=3)
  Y=matrix(rnorm(1000*30),nrow=1000,ncol=30)
  colnames(X)=paste0("X",1:ncol(X))
  colnames(Y)=paste0("Y",1:ncol(Y))

  at=association_testing_quantile(Y, X, Z)
  
  expect_false(any(is.na(at)))
  
  at=association_testing_quantile(Y, X, Z, K=1)
  
  expect_false(any(is.na(at)))

})
