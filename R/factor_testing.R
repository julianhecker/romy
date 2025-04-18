



#' ROMY factor association testing
#'
#' This function performs factor association testing between variables in Y and an element of X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param num_factors Number of factors to test, default is 3.
#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param prediction_x Prediction model to predict X given Z. Default is linear regression.
#' @param learner_y Prediction model function to learn prediction of Y given Z. Default is linear regression.
#' @param prediction_y Prediction model to predict Y given Z. Default is linear regression.
#' @param method Method for controlling cross-fitting. Options are 'double_cf', 'single_cf', or 'no_split'. Default is 'double_cf'.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param split_ratio Ratio for splitting of training data for the two prediction tasks. Default is 0.5:0.5.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation. Default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
association_factor_testing=function(Y, X, Z, num_factors=3, 
learner_x=lm_learner, prediction_x=lm_predict, learner_y=lm_learner, prediction_y=lm_predict, method="double_cf",
K=5, split_ratio=c(0.3333, 0.3333, 1-0.3333-0.3333), parallel=FALSE, BPPARAM=NULL)
{
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, nomissX=FALSE, nomissY=FALSE, parallel=parallel, BPPARAM=BPPARAM)
	  
	  .check_model(learner=learner_x, prediction=prediction_x)
	  .check_model(learner=learner_y, prediction=prediction_y)
	  
	  if(!(method %in% c("double_cf","single_cf","no_split")))
	  {
		 stop("method unknown.")
	  }
	  
	  
	  ################################################
	  ### get dimensions
	  N=nrow(Y)
	  n_X=ncol(X)
	  n_Y=ncol(Y)
	  
	  
	  
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
	  if(method=="double_cf")
	  {
		splits=.create_splits_3subs(K=K, N=N, split_ratio=split_ratio)
	  }
	  if(method=="single_cf")
	  {
		splits=.create_splits_3subs(K=K, N=N, split_ratio=split_ratio, single=TRUE)
	  }
	  if(method=="no_split")
	  {
		splits=.create_no_split_data(N=N, subs=3)
	  }
	  ################################################################
	  ### factor identification
	  factors=.factor_estimation(Y=Y, X=X, Z=Z, splits=splits, training_part=1, num_factors=num_factors)
	  
	  ################################################################
	  ### compute X residuals
	  if(!parallel){
		  Xrs=.get_residuals(Outcome=X, Z=Z, splits=splits, training_part=2, learner=learner_x, prediction=prediction_x)
	  }
	  if(parallel){
		  Xrs=.get_residuals_parallel(Outcome=X, Z=Z, splits=splits, BPPARAM=BPPARAM, training_part=2, learner=learner_x, prediction=prediction_x)
	  }
	  Xrs=Xrs[[1]]
	  
	    
	  ###################################################################
	  ### perform double machine learning testing
	  
	  results_df=.dml_association_factor_testing(Xrs=Xrs, Y=Y, Z=Z, col_X=colnames(X), splits=splits, factors, num_factors = 3)
	  
	  
	  ###################################################################
	  ### 
	  colnames(results_df)[1] <- c("X_id")
	  
	  if(ncol(results_df)>3){
	    results_df[, 2:ncol(results_df)] <- lapply(results_df[, 2:ncol(results_df)], as.numeric)
	    
	    ##if m==1 => p_overall and p_ACAT are not necessary
	    results_df$p_ACAT <- apply(results_df[, 2:ncol(results_df)], 1, function(x){if(any(is.na(x))){return(NA)}else{return(ACAT(x))}})
	  }
	  if(ncol(results_df)==3){
	    results_df[,2]=as.numeric(results_df[,2])
	    results_df=results_df[,1:2]
	    colnames(results_df)[2]='p_value'
	  }
	  
	  return(results_df)
}

#' This function performs the factor estimation
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param splits The splits of the samples determining the K-fold cross fitting.#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param num_factors Number of factors to test, default is 3.
#'
.factor_estimation=function(Y, X, Z, splits, training_part=1, num_factors=3){
  
  K=length(splits)
  factors=list()
  for(p in 1:ncol(X))
  {
      sub_factors=list()
      for(k in 1:K)
      {
        inds_train=splits[[k]]$inds_train[[training_part]]
        Y_train=Y[inds_train,]
        Y_train_p=Y[inds_train,]
        predictors=as.matrix(cbind(X[inds_train,p], Z[inds_train,]))
        XXi = solve(t(predictors) %*% predictors)
        for(i in 1:ncol(Y))
        {
           betas=XXi %*% t(predictors) %*% Y_train[,i]
           Y_train_p[,i]=X[inds_train,p] * betas[1]
        }
        evd=eigen(t(Y_train_p) %*% Y_train_p)
        v=evd$vectors[,1:num_factors]
        for(n in 1:num_factors)
        {
          if(v[1,n]<0){v[,n]=(-1)*v[,n]}
        }
        sub_factors[[k]]=v
      }
      factors[[p]]=sub_factors
  }
  return(factors)
}

#' Covariate adjustment and residual computation
#'
#' This function performs the covariate adjustment for the factor association testing
#'
#' @param Y A matrix containing the Y measurements.
#' @param Z A matrix containing the covariates Z.
#' @param ind_x Index of element of X to be tested.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param factors Factors estimated from data.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Function pointer to prediction model (prediction). Default is 'lm_predict_simple'.
#' @param num_factors Number of factors to test, default is 3.
#'
.get_individual_residuals=function(Y, Z, ind_x, splits, factors, training_part=3, learner=lm_learner_simple, prediction=lm_predict_simple, num_factors=3)
{
  
  K=length(splits)
  Residuals=matrix(0, nrow=nrow(Y), ncol=num_factors)
  
  for(k in 1:K)
  {
        inds_train=splits[[k]]$inds_train[[training_part]] 
        inds_test=splits[[k]]$inds_test
        
        O_train=Y[inds_train,];
        Z_train=Z[inds_train,];
        
        O_test=Y[inds_test,];
        Z_test=Z[inds_test,];
        
        tmp_v=factors[[ind_x]][[k]]
        
        O_train = O_train %*% tmp_v
        O_test = O_test %*% tmp_v
        
        data_Z_train=data.frame(Z=Z_train)
        colnames(data_Z_train)=paste0("Z", 1:ncol(Z_train))
        
        data_Z_test=data.frame(Z=Z_test)
        colnames(data_Z_test)=paste0("Z", 1:ncol(Z_test))
        
        for(j in 1:num_factors)
        {
          O_train_j = O_train[,j]
          O_test_j = O_test[,j]
          fit=learner(O_train_j, data_Z_train)
          res=O_test_j-prediction(fit, data_Z_test)  
          Residuals[inds_test, j]=res
        }
        
        
        
        
  }
   
  return(Residuals)
}

#' This function computes the test statistics for double machine learning 
#' based association testing, given the residual measurements and index pairs.
#' using parallel computations
#'
#' @param Xrs Residual matrices for X.
#' @param Y Matrix of Y measurements.
#' @param Z Matrix of Z measurements.
#' @param col_X Column names for X.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param factors Factors estimated from data.
#' @param num_factors Number of factors to test, default is 3.
#'
.dml_association_factor_testing=function(Xrs, Y, Z, col_X, splits, factors, num_factors=3){

  indices_x=as.matrix(1:length(col_X))
  K=length(splits)
	results <- apply(indices_x, 1, function(idx){
	  
      		i <- idx
      		stats <- numeric(num_factors)
      		variances <- numeric(num_factors)
      		covariance_mat <- matrix(0, nrow = num_factors, ncol = num_factors)
      
      		Yrs=.get_individual_residuals(Y=Y, Z=Z, ind_x=i, splits=splits, factors=factors, training_part=3, num_factors=num_factors)
      		
      		stats=colSums(Xrs[,i]*Yrs)
			stats_means=colMeans(Xrs[,i]*Yrs)
			tmp=Xrs[,i]*Yrs-stats_means
      	
      		covariance_mat = t(tmp) %*% (tmp)
			variances=diag(covariance_mat)
      		
      		p <- 2 * pnorm(-abs(stats / sqrt(variances)), 0, 1) # compute asymptotic p-values based on normal distribution
      		
      		p_overall=NA
      		tryCatch({
      			stat_overall=t(stats) %*% solve(covariance_mat) %*% stats
      		    p_overall=pchisq(stat_overall, df=num_factors, lower.tail=FALSE) # compute overall asymptotic p-value based on chisq distribution
      			
      		}, error=function(e){
      				warning("Covariance matrix not invertible, returning NA.")
      				
      		})
      		
      
      		c(matrix1_col = col_X[i], p_values=p, p_overall=p_overall)
	  })
	  results<- as.data.frame(t(results))
	  return(results)
}


