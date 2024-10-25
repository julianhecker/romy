#' ROMY component testing
#'
#'
#' This function performs component testing between Y and components of X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A list of matrices, containing the measurements for X as separate multi-dimensional components.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param learner_x Prediction model function to learn prediction of X given Z. Default is linear regression.
#' @param prediction_x Prediction model to predict X given Z. Default is linear regression.
#' @param learner_y Prediction model function to learn prediction of Y given Z. Default is linear regression.
#' @param prediction_y Prediction model to predict Y given Z. Default is linear regression.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param split_ratio Ratio for splitting of training data for the two prediction tasks. Default is 0.5:0.5.
#' @param num_factors Number of factors to be extracted in Singular Value Decomposition. Default is 3.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation. Default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
component_testing=function(Y, X, Z, learner_x=lm_learner, prediction_x=lm_predict, 
learner_y=lm_learner, prediction_y=lm_predict, K=5, split_ratio=c(0.3333,0.3333,1-2*0.3333), num_factors=3, parallel=FALSE, BPPARAM=NULL)
{
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, split_ratio=split_ratio, num_splits=3, nomissX=TRUE, nomissY=FALSE, K=K, parallel=parallel, BPPARAM=BPPARAM)
	  .check_model(learner=learner_x, prediction=prediction_x)
	  .check_model(learner=learner_y, prediction=prediction_y)
	  ################################################
	  ### get dimensions
	  
	  N=nrow(Y)
	  ### assume Y has 1 column
	  ################################################################
	  ### scale data, missing data? 
	  Y=scale(Y, center = TRUE, scale = FALSE) # ignores NA value
	  X=scale(X, center = TRUE, scale = TRUE)
	  Z=scale(Z, center = TRUE, scale = TRUE)
	  if(is.null(colnames(X))){
	     colnames(X)=paste0("X",1:ncol(X))
	  }
	  if(is.null(colnames(Y))){
	     colnames(Y)=paste0("Y",1:ncol(Y))
	  }
	  
	  ################################################################
	  ### get data splits
	  
	  splits=.create_splits_3subs(K=K, N=N, split_ratio=split_ratio)
	  
	  ################################################################
	  ### perform factor analysis 
	  
	  F_p=.factor_analysis(X=X, splits=splits, training_part=1, num_factors=num_factors)

	  
	  ################################################################
	  ### perform component testing
	  if(!parallel){
		  Y_resid=.get_residuals(Outcome=Y, Z=Z, splits=splits, training_part=2, learner=learner_x, prediction=prediction_x)
		  F_resid=.get_residuals(Outcome=F_p, Z=Z, splits=splits, training_part=3, learner=learner_x, prediction=prediction_x)
	  }
	  if(parallel){
		  Y_resid=.get_residuals_parallel(Outcome=Y, Z=Z, splits=splits, BPPARAM=BPPARAM, training_part=2, learner=learner_x, prediction=prediction_x)
		  F_resid=.get_residuals_parallel(Outcome=F_p, Z=Z, splits=splits, BPPARAM=BPPARAM, training_part=3, learner=learner_x, prediction=prediction_x)
	  }
	  pval=.dml_association_testing_comp(Yrs=Y_resid, Xrs=F_resid)
	  
	  ###################################################################
	  return(list(pval=pval, F_p=F_p))
}


#' This function computes the test statistic for double machine learning based association testing, given the residual measurements.
#'
#' @param Xrs Residual matrices for X.
#' @param Yrs Residual matrices for Y.
#'
.dml_association_testing_comp=function(Yrs, Xrs){

	m1=ncol(Xrs[[1]])
	m2=ncol(Yrs[[1]])
	m=m1*m2
		
	stats <- numeric(m)
	covariance_mat <- matrix(0, nrow = m, ncol = m)

	ctr1 <- 0

	for (k1 in 1:m1) {
		X_k1 <- Xrs[[1]][, k1]  
		  
		for (l1 in 1:m2) {
			ctr1 <- ctr1 + 1
			Y_l1 <- Yrs[[1]][, l1]  
						
						# Calculate stat and variances
			stats[ctr1] <- sum(X_k1 * Y_l1, na.rm = TRUE)
			
						
			ctr2 <- 0
						
			for (k2 in 1:m1) {
				X_k2 <- Xrs[[1]][, k2] 
						  
				for (l2 in 1:m2) {
					ctr2 <- ctr2 + 1
					Y_l2 <- Yrs[[1]][, l2]  
								
					covariance_mat[ctr1, ctr2] <- sum(X_k1 * Y_l1 * X_k2 * Y_l2, na.rm = TRUE)
				}
			}
		}
	}
	tryCatch({
			stat_overall=t(stats) %*% solve(covariance_mat) %*% stats
		    p_overall=pchisq(stat_overall, df=m, lower.tail=FALSE) # compute overall asymptotic p-value based on chisq distribution
			
	}, error=function(e){
				warning("Covariance matrix not invertible, returning NA.")
				p_overall=NA
	})
	

	return(as.numeric(p_overall))
}



.factor_analysis=function(X, splits, training_part,num_factors){

    K=length(splits)
	N=nrow(X)
	F_p=matrix(0, nrow=N, ncol=num_factors)
	
	for(k in 1:K){
	
			inds_train=splits[[k]]$inds_train[[training_part]] 
			inds_test=splits[[k]]$inds_test
			
			
			X_train=as.matrix(X[inds_train,]);
			
			X_test=as.matrix(X[inds_test,]);
			
			
			svd_X=svd(X_train)

			v=svd_X$v[,1:num_factors]
			for(i in 1:num_factors)
			{
			  v[,i]=v[,i]/svd_X$d[i] * sqrt(N)
			}

			tmp_F=X_test %*% v

			F_p[inds_test,]=tmp_F

	}
	return(F_p)
}