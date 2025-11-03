library(romy)

test_that("check_inter", {
	X=matrix(rnorm(1000*30),nrow=1000,ncol=30)
	Z=matrix(rnorm(1000*5),nrow=1000,ncol=5)
	Y=matrix(rnorm(1000*1),nrow=1000,ncol=1)
	colnames(X)=paste0("X",1:ncol(X))
	colnames(Y)=paste0("Y",1:ncol(Y))
	colnames(Z)=paste0("Z",1:ncol(Z))

	inter=romy_inter(Y=Y, X=X, Z=Z)
  
    expect_false(any(is.na(inter)))

})