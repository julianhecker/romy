#' ROMY association testing
#'
#' This function performs association testing between variables in Y and X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param index_pairs A dataframe with two columns containing the pairs of indices in X and Y to be tested, respectively. Default is NULL, leading to testing all possible pairs.
#' @param transform_functions_x List of functions to transform the measurements for X. All transformations will be considered in testing.
#' @param transform_functions_y List of functions to transform the measurements for Y. All transformations will be considered in testing.
#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param prediction_x Prediction model to predict X given Z. Default is linear regression.
#' @param learner_y Prediction model function to learn prediction of Y given Z. Default is linear regression.
#' @param prediction_y Prediction model to predict Y given Z. Default is linear regression.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param split_ratio Ratio for splitting of training data for the two prediction tasks. Default is 0.5:0.5.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation. Default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
association_testing=function(Y, X, Z, index_pairs=NULL, 
transform_functions_x=NULL, transform_functions_y=NULL, 
learner_x=lm_learner, prediction_x=lm_predict, learner_y=lm_learner, prediction_y=lm_predict, 
K=5, split_ratio=c(0.5,0.5), parallel=FALSE, BPPARAM=NULL)
{
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, split_ratio=split_ratio, num_splits=2, nomissX=FALSE, nomissY=FALSE, K=K, parallel=parallel, BPPARAM=BPPARAM)
	  
	  .check_model(learner=learner_x, prediction=prediction_x)
	  .check_model(learner=learner_y, prediction=prediction_y)
	  
	  if(is.null(transform_functions_x)){
		 fun_x <- function(x) { return(x) }
		 transform_functions_x=list()
		 transform_functions_x[[1]]=fun_x
	  }
	  if(is.null(transform_functions_y)){
		 fun_y <- function(y) { return(y) }
		 transform_functions_y=list()
		 transform_functions_y[[1]]=fun_y
      }
	  
	  ################################################
	  ### get dimensions
	  N=nrow(Y)
	  n_X=ncol(X)
	  n_Y=ncol(Y)
	  if(!is.null(index_pairs) & (!all(1:n_Y %in% unique(index_pairs[,2])) | (!all(1:n_X %in% unique(index_pairs[,1]))))){stop("index_pairs contains not all components of Y and X")}
	  if(is.null(index_pairs)){ index_pairs <- expand.grid(1:n_X, 1:n_Y)}
	  
	  ################################################################
	  ### scale data, missing data?
	  
	  Y=scale(Y, center = TRUE, scale = FALSE) # ignores NA values
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
	  splits=.create_splits_2subs(K=K, N=N, split_ratio=split_ratio)
	  
	  ################################################################
	  ### compute residuals
	  if(!parallel){
		  Xrs=.get_residuals(Outcome=X, Z=Z, transform_functions=transform_functions_x, splits=splits, training_part=1, learner=learner_x, prediction=prediction_x)
		  Yrs=.get_residuals(Outcome=Y, Z=Z, transform_functions=transform_functions_y, splits=splits, training_part=2, learner=learner_y, prediction=prediction_y)
	  }
	  if(parallel){
		  Xrs=.get_residuals_parallel(Outcome=X, Z=Z, transform_functions=transform_functions_x, splits=splits, BPPARAM=BPPARAM, training_part=1, learner=learner_x, prediction=prediction_x)
		  Yrs=.get_residuals_parallel(Outcome=Y, Z=Z, transform_functions=transform_functions_y, splits=splits, BPPARAM=BPPARAM, training_part=2, learner=learner_y, prediction=prediction_y)
	  }
	  ###################################################################
	  ### perform double machine learning testing
	  if(!parallel){
		results_df=.dml_association_testing(Xrs=Xrs, Yrs=Yrs, col_X=colnames(X), col_Y=colnames(Y), index_pairs=index_pairs)
	  }
	  if(parallel){
		results_df=.dml_association_testing_parallel(Xrs=Xrs, Yrs=Yrs, col_X=colnames(X), col_Y=colnames(Y), index_pairs=index_pairs, BPPARAM=BPPARAM)
	  }
	  ###################################################################
	  ### add ACAT if multiple tests performed, formatting
	  
	  colnames(results_df)[1:2] <- c("X_id", "Y_id")
	  
	  if(ncol(results_df)>4){
		  results_df[, 3:ncol(results_df)] <- lapply(results_df[, 3:ncol(results_df)], as.numeric)
		
		  ##if m==1 => p_overall and p_ACAT are not necessary
		  results_df$p_ACAT <- apply(results_df[, 3:ncol(results_df)], 1, function(x){if(any(is.na(x))){return(NA)}else{return(ACAT(x))}})
	  }
	  if(ncol(results_df)==4){
		  results_df[,3]=as.numeric(results_df[,3])
		  results_df=results_df[,1:3]
		  colnames(results_df)[3]='p_value'
	  }
	  return(results_df)
}





#' This function computes the test statistics for double machine learning 
#' based association testing, given the residual measurements and index pairs.
#'
#' @param Xrs Residual matrices for X.
#' @param Yrs Residual matrices for Y.
#' @param col_X Column names for X.
#' @param col_Y Column names for Y.
#' @param index_pairs Pairs of X-Y to test.
#'
.dml_association_testing=function(Xrs, Yrs, col_X, col_Y, index_pairs){

	results <- apply(index_pairs, 1, function(idx){
		i <- idx[1]
		j <- idx[2]
		m1=length(Xrs)
		m2=length(Yrs)
		m=m1*m2
		
		stats <- numeric(m)
		variances <- numeric(m)
		covariance_mat <- matrix(0, nrow = m, ncol = m)

		ctr1 <- 0

		for (k1 in 1:m1) {
			  X_k1 <- Xrs[[k1]][, i]  
		  
			  for (l1 in 1:m2) {
					ctr1 <- ctr1 + 1
					Y_l1 <- Yrs[[l1]][, j]  
					
					# Calculate stat and variances
					stats[ctr1] <- sum(X_k1 * Y_l1, na.rm = TRUE)
					variances[ctr1] <- sum(X_k1^2 * Y_l1^2, na.rm = TRUE)
					
					ctr2 <- 0
					
					for (k2 in 1:m1) {
						  X_k2 <- Xrs[[k2]][, i] 
					  
						  for (l2 in 1:m2) {
							ctr2 <- ctr2 + 1
							Y_l2 <- Yrs[[l2]][, j]  
							
							covariance_mat[ctr1, ctr2] <- sum(X_k1 * Y_l1 * X_k2 * Y_l2, na.rm = TRUE)
						  }
					}
			  }
		}

		p <- 2 * pnorm(-abs(stats / sqrt(variances)), 0, 1) # compute asymptotic p-values based on normal distribution
		
		p_overall=NA
		tryCatch({
			stat_overall=t(stats) %*% solve(covariance_mat) %*% stats
		    p_overall=pchisq(stat_overall, df=m, lower.tail=FALSE) # compute overall asymptotic p-value based on chisq distribution
			
		}, error=function(e){
				warning("Covariance matrix not invertible, returning NA.")
				
		})
		

		c(matrix1_col = col_X[i], matrix2_col = col_Y[j], p_values=p, p_overall=p_overall)
	  })
	  results<- as.data.frame(t(results))
	  return(results)
}


#' This function computes the test statistics for double machine learning 
#' based association testing, given the residual measurements and index pairs.
#' using parallel computations
#'
#' @param Xrs Residual matrices for X.
#' @param Yrs Residual matrices for Y.
#' @param col_X Column names for X.
#' @param col_Y Column names for Y.
#' @param index_pairs Pairs of X-Y to test.
#' @param BPPARAM BPPARAM object for parallel computing.
#'
.dml_association_testing_parallel=function(Xrs, Yrs, col_X, col_Y, index_pairs, BPPARAM){

	compute_results <- function(id){
	
		i <- index_pairs[id,1]
		j <- index_pairs[id,2]
		m1=length(Xrs)
		m2=length(Yrs)
		m=m1*m2
		
		stats <- numeric(m)
		variances <- numeric(m)
		covariance_mat <- matrix(0, nrow = m, ncol = m)

		ctr1 <- 0

		for (k1 in 1:m1) {
			  X_k1 <- Xrs[[k1]][, i]  
		  
			  for (l1 in 1:m2) {
					ctr1 <- ctr1 + 1
					Y_l1 <- Yrs[[l1]][, j]  
					
					# Calculate stat and variances
					stats[ctr1] <- sum(X_k1 * Y_l1, na.rm = TRUE)
					variances[ctr1] <- sum(X_k1^2 * Y_l1^2, na.rm = TRUE)
					
					ctr2 <- 0
					
					for (k2 in 1:m1) {
						  X_k2 <- Xrs[[k2]][, i] 
					  
						  for (l2 in 1:m2) {
							ctr2 <- ctr2 + 1
							Y_l2 <- Yrs[[l2]][, j]  
							
							covariance_mat[ctr1, ctr2] <- sum(X_k1 * Y_l1 * X_k2 * Y_l2, na.rm = TRUE)
						  }
					}
			  }
		}

		p <- 2 * pnorm(-abs(stats / sqrt(variances)), 0, 1) # compute asymptotic p-values based on normal distribution
		
		p_overall=NA
		tryCatch({
			stat_overall=t(stats) %*% solve(covariance_mat) %*% stats
		    p_overall=pchisq(stat_overall, df=m, lower.tail=FALSE) # compute overall asymptotic p-value based on chisq distribution
			
		}, error=function(e){
				warning("Covariance matrix not invertible, returning NA.")
				
		})
		

		c(matrix1_col = col_X[i], matrix2_col = col_Y[j], p_values=p, p_overall=p_overall)
	  }
	  
	  #results_list <- BiocParallel::bplapply(seq_len(nrow(index_pairs)), compute_results, BPPARAM = BPPARAM)
	  num_cores = BiocParallel::multicoreWorkers()
	  results = BiocParallel::bplapply(
			split(seq_len(nrow(index_pairs)), cut(seq_len(nrow(index_pairs)), breaks=num_cores)), 
			function(ids) lapply(ids, compute_results), 
			BPPARAM = BPPARAM
	  )
	  results = do.call(c, results)
	  
	  results=do.call(rbind, results)
	  results=as.data.frame(results)
	  return(results)
}