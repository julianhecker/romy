library(romy)

test_that("check_component_testing", {
	X=matrix(rnorm(1000*30),nrow=1000,ncol=30)
	Z=matrix(rnorm(1000*5),nrow=1000,ncol=5)
	Y=matrix(rnorm(1000),nrow=1000,ncol=1)
	colnames(Y)=paste0("Y",1:ncol(Y))
	colnames(Z)=paste0("Z",1:ncol(Z))

    ct=component_testing(Y, X, Z)
  
    expect_false(any(is.na(ct$pval)))

})