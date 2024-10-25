library(romy)

test_that("check_covar_testing", {
    X=matrix(rnorm(1000*3),nrow=1000,ncol=3)
	Z=matrix(rnorm(1000*3),nrow=1000,ncol=3)
	Y=matrix(rnorm(1000*10),nrow=1000,ncol=10)
	colnames(X)=paste0("X",1:ncol(X))
	colnames(Y)=paste0("Y",1:ncol(Y))

	ct=covar_testing(Y, X, Z)
  
    expect_false(any(is.na(ct)))

})
