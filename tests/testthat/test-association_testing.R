library(romy)

test_that("check_association_testing", {
  X=matrix(rnorm(1000*3),nrow=1000,ncol=3)
  Z=matrix(rnorm(1000*3),nrow=1000,ncol=3)
  Y=matrix(rnorm(1000*30),nrow=1000,ncol=30)
  colnames(X)=paste0("X",1:ncol(X))
  colnames(Y)=paste0("Y",1:ncol(Y))

  at=association_testing(Y, X, Z)
  
  expect_false(any(is.na(at)))
  
  at=association_testing(Y, X, Z, method='no_split', K=1)
  
  expect_false(any(is.na(at)))

})
