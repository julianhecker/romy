#' ROMY-CIT
#'
#' This function performs association testing between variables in Y and X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z, for adjusting Z. No missing values allowed.
#' @param Z2 A matrix containing the measurements for Z2, for adjusting X. No missing values allowed, default is NULL. In this case, the computations use Z2=Z.
#' @param index_pairs A dataframe with two columns containing the pairs of indices in X and Y to be tested, respectively. Default is NULL, leading to testing all possible pairs.
#' @param transform_functions_x List of functions to transform the measurements for X. All transformations b(x) will be considered in testing. Default is NULL, implying b(x)=x.
#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param prediction_x Prediction model to predict X given Z. Default is linear regression.
#' @param learner_y Prediction model function to learn prediction of Y given Z. Default is linear regression.
#' @param prediction_y Prediction model to predict Y given Z. Default is linear regression.
#' @param method Method for controlling cross-fitting. Options are 'double_cf', 'single_cf', or 'no_split'. Default is 'no_split'.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation. Default is NULL.
#'
#' return Dataframe containing the p-values for all tests, including ACAT results (if more than 1 test per combination).
#'
#' @export
romy_cit=function(Y, X, Z, Z2=NULL, index_pairs=NULL, 
transform_functions_x=NULL, 
learner_x=lm_learner, prediction_x=lm_predict, learner_y=lm_learner, prediction_y=lm_predict, method="no_split",
K=5, parallel=FALSE, BPPARAM=NULL)
{
	  if(is.null(Z2))
	  {
	      Z2=Z
	  }
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, Z2=Z2, nomissX=FALSE, nomissY=FALSE, parallel=parallel, BPPARAM=BPPARAM)
	  
	  .check_model(learner=learner_x, prediction=prediction_x)
	  .check_model(learner=learner_y, prediction=prediction_y)
	  
	  if(!(method %in% c("double_cf","single_cf","no_split")))
	  {
	  	 stop("method not implemented for association testing.")
	  }
	  
	  if(is.null(transform_functions_x)){
		 fun_x <- function(x) { return(x) }
		 transform_functions_x=list()
		 transform_functions_x[[1]]=fun_x
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
	  Z2=scale(Z2, center = TRUE, scale = TRUE)
	  
	  if(is.null(colnames(X))){
	     colnames(X)=paste0("X",1:ncol(X))
	  }
	  if(is.null(colnames(Y))){
	     colnames(Y)=paste0("Y",1:ncol(Y))
	  }
	  ################################################################
	  ### get data splits
	  if(method=="double_cf")
	  {
		splits=.create_splits(N=N, K=K, subs=2)
		training_xrs=1
		training_yrs=2
	  }
	  if(method=="single_cf")
	  {
		splits=.create_splits(N=N, K=K, subs=1)
		training_xrs=1
		training_yrs=1
	  }
	  if(method=="no_split")
	  {
		# K ignored and set to K=1
		splits=.create_splits(N=N, K=1, subs=1)
		training_xrs=1
		training_yrs=1
	  }

	  
	  ################################################################
	  ### compute residuals
	  if(!parallel){
		  Xrs=.get_residuals(Outcome=X, Z=Z2, transform_functions=transform_functions_x, splits=splits, training_part=training_xrs, learner=learner_x, prediction=prediction_x)
		  Yrs=.get_residuals(Outcome=Y, Z=Z, splits=splits, training_part=training_yrs, learner=learner_y, prediction=prediction_y)
	  }
	  if(parallel){
		  Xrs=.get_residuals_parallel(Outcome=X, Z=Z2, transform_functions=transform_functions_x, splits=splits, BPPARAM=BPPARAM, training_part=training_xrs, learner=learner_x, prediction=prediction_x)
		  Yrs=.get_residuals_parallel(Outcome=Y, Z=Z, splits=splits, BPPARAM=BPPARAM, training_part=training_yrs, learner=learner_y, prediction=prediction_y)
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
		m=length(Xrs)
		
		
		stats <- numeric(m)
		variances <- numeric(m)
		covariance_mat <- matrix(0, nrow = m, ncol = m)
		
		N=nrow(Yrs[[1]])
		
		tmp_mat=matrix(0, nrow=N, ncol=m)

		for (k in 1:m) {
			 X_ki <- Xrs[[k]][, i]  
		  
			  
			 Y_j <- Yrs[[1]][, j]  
					
			 tmp_mat[,k]=X_ki * Y_j
					
					
			  
		}
		stats=colSums(tmp_mat, na.rm=T)
		tmp_mat=tmp_mat-colMeans(tmp_mat, na.rm=T)
		variances=colSums(tmp_mat**2, na.rm=T)
		covariance_mat=.na_safe_matrix_product(t(tmp_mat), tmp_mat)

		p <- exp(log(2) + pnorm(-abs(stats / sqrt(variances)), 0, 1, log.p = TRUE))
		
		p_overall=NA
		tryCatch({
			stat_overall=t(stats) %*% solve(covariance_mat) %*% stats
			p_overall <- exp(pchisq(stat_overall, df = m, lower.tail = FALSE, log.p = TRUE))
			
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
		m=length(Xrs)
		
		
		stats <- numeric(m)
		variances <- numeric(m)
		covariance_mat <- matrix(0, nrow = m, ncol = m)

		N=nrow(Yrs[[1]])
		
		tmp_mat=matrix(0, nrow=N, ncol=m)

	
		for (k in 1:m) {
			  X_ki <- Xrs[[k]][, i]  
		  
			  
			  Y_j <- Yrs[[1]][, j]  
					
			  tmp_mat[,k]=X_ki * Y_j
					
					
			  
		}
		stats=colSums(tmp_mat, na.rm=T)
		tmp_mat=tmp_mat-colMeans(tmp_mat, na.rm=T)
		variances=colSums(tmp_mat**2, na.rm=T)
		covariance_mat=.na_safe_matrix_product(t(tmp_mat), tmp_mat)

		p <- 2 * pnorm(-abs(stats / sqrt(variances)), 0, 1) # compute asymptotic p-values based on normal distribution
		
		p_overall=NA
		tryCatch({
			stat_overall=t(stats) %*% solve(covariance_mat) %*% stats
			p_overall <- exp(pchisq(stat_overall, df = m, lower.tail = FALSE, log.p = TRUE))
			
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