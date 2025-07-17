



#' ROMY-FCIT
#'
#' This function performs factor association testing between variables in Y and an element of X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z, for adjusting Z. No missing values allowed.
#' @param Z2 A matrix containing the measurements for Z2, for adjusting X. No missing values allowed, default is NULL. In this case, Z2=Z.
#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param prediction_x Prediction model to predict X given Z. Default is linear regression.
#' @param learner_y Prediction model function to learn prediction of Y given Z. Default is linear regression.
#' @param prediction_y Prediction model to predict Y given Z. Default is linear regression.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation. Default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
romy_fcit=function(Y, X, Z, Z2=NULL,
learner_x=lm_learner, prediction_x=lm_predict, learner_y=lm_learner, prediction_y=lm_predict, 
K=5, parallel=FALSE, BPPARAM=NULL)
{
	  if(is.null(Z2))
	  {
	      Z2=Z
	  }
	  ### initial checks
	  #.input_checks(Y=Y, X=X, Z=Z, Z2=Z2, nomissX=TRUE, nomissY=TRUE, parallel=parallel, BPPARAM=BPPARAM)
	  
	  #.check_model(learner=learner_x, prediction=prediction_x)
	  #.check_model(learner=learner_y, prediction=prediction_y)
	  
	  
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
	  Z2=scale(Z2, center = TRUE, scale = TRUE)
	  
	  if(is.null(colnames(X))){
	     colnames(X)=paste0("X",1:ncol(X))
	  }
	  if(is.null(colnames(Y))){
	     colnames(Y)=paste0("Y",1:ncol(Y))
	  }
	  ################################################################
	  ### get data splits
	  
	  
	  splits=.create_splits(N=N, K=K, subs=3)
	  training_factor=1
	  training_x=2
	  training_y=3
	  
	  ################################################################
	  ### factor identification
	  if(!parallel){
		factors=.factor_estimation(Y=Y, X=X, Z=Z, splits=splits, training_part=training_factor)
	  }
	  if(parallel){
		factors=.factor_estimation_parallel(Y=Y, X=X, Z=Z, splits=splits, training_part=training_factor, BPPARAM=BPPARAM)
	  }
	  ################################################################
	  ### compute X residuals
	  if(!parallel){
		  Xrs=.get_residuals(Outcome=X, Z=Z2, splits=splits, training_part=training_x, learner=learner_x, prediction=prediction_x)
	  }
	  if(parallel){
		  Xrs=.get_residuals_parallel(Outcome=X, Z=Z2, splits=splits, BPPARAM=BPPARAM, training_part=training_x, learner=learner_x, prediction=prediction_x)
	  }
	  Xrs=Xrs[[1]]
	  
	    
	  ###################################################################
	  ### perform double machine learning testing
	  if(!parallel){
		results_df=.dml_association_factor_testing(Xrs=Xrs, Y=Y, Z=Z, col_X=colnames(X), splits=splits, factors=factors, training_y=training_y)
	  }
	  if(parallel){
		results_df=.dml_association_factor_testing_parallel(Xrs=Xrs, Y=Y, Z=Z, col_X=colnames(X), splits=splits, factors=factors, BPPARAM=BPPARAM,
		training_y=training_y)
	  }
	  
	  ###################################################################
	  ### 
	  colnames(results_df)[1] <- c("X_id")
	  
	  
	  results_df[,2]=as.numeric(results_df[,2])
	  results_df=results_df[,1:2]
	  colnames(results_df)[2]='p_value'
	  
	  
	  return(results_df)
}

#' This function performs the factor estimation
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param splits The splits of the samples determining the K-fold cross fitting.#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#'
.factor_estimation=function(Y, X, Z, splits, training_part=1){
  
  K=length(splits)
  factors=list()
  for(p in 1:ncol(X))
  {
      sub_factors=list()
      for(k in 1:K)
      {
        inds_train=splits[[k]]$inds_train[[training_part]]
        Y_train=as.matrix(Y[inds_train,])
        Y_train_p=as.matrix(Y[inds_train,])
        predictors=as.matrix(cbind(X[inds_train,p], Z[inds_train,]))
        XXi = solve(t(predictors) %*% predictors)
        for(i in 1:ncol(Y))
        {
           betas=XXi %*% t(predictors) %*% Y_train[,i]
           Y_train_p[,i]=X[inds_train,p] * betas[1]
        }
        evd=eigen(t(Y_train_p) %*% Y_train_p)
        v=as.matrix(evd$vectors[,1])
        
        if(v[1,1]<0){v[,1]=(-1)*v[,1]}
        
        sub_factors[[k]]=v
      }
      factors[[p]]=sub_factors
  }
  return(factors)
}


#' This function performs the factor estimation
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param splits The splits of the samples determining the K-fold cross fitting.#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation.
#'
.factor_estimation_parallel=function(Y, X, Z, splits, training_part=1, BPPARAM){
  
  K=length(splits)
  factors=list()
  compute_results=function(id){
      sub_factors=list()
      for(k in 1:K)
      {
        inds_train=splits[[k]]$inds_train[[training_part]]
        Y_train=as.matrix(Y[inds_train,])
        Y_train_p=as.matrix(Y[inds_train,])
        predictors=as.matrix(cbind(X[inds_train,id], Z[inds_train,]))
        XXi = solve(t(predictors) %*% predictors)
        for(i in 1:ncol(Y))
        {
           betas=XXi %*% t(predictors) %*% Y_train[,i]
           Y_train_p[,i]=X[inds_train,id] * betas[1]
        }
        evd=eigen(t(Y_train_p) %*% Y_train_p)
        v=as.matrix(evd$vectors[,1])
        
        if(v[1,1]<0){v[,1]=(-1)*v[,1]}
        
        sub_factors[[k]]=v
      }
      c(sub_factors)
  }
  
  num_cores = BiocParallel::multicoreWorkers()
  factors = BiocParallel::bplapply(
			split(seq_len(ncol(X)), cut(seq_len(ncol(X)), breaks=num_cores)), 
			function(ids) lapply(ids, compute_results), 
			BPPARAM = BPPARAM
  )
  factors = do.call(c, factors)
	
  
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
#'
.get_individual_residuals=function(Y, Z, ind_x, splits, factors, training_part=3, learner=lm_learner_simple, prediction=lm_predict_simple)
{
  
  K=length(splits)
  Residuals=matrix(0, nrow=nrow(Y), ncol=1)
  
  for(k in 1:K)
  {
        inds_train=splits[[k]]$inds_train[[training_part]] 
        inds_test=splits[[k]]$inds_test
        
        O_train=as.matrix(Y[inds_train,]);
        Z_train=as.matrix(Z[inds_train,]);
        
        O_test=as.matrix(Y[inds_test,]);
        Z_test=as.matrix(Z[inds_test,]);
        
        tmp_v=as.matrix(factors[[ind_x]][[k]])
        
        O_train = O_train %*% tmp_v
        O_test = O_test %*% tmp_v
        
        data_Z_train=data.frame(Z=Z_train)
        colnames(data_Z_train)=paste0("Z", 1:ncol(Z_train))
        
        data_Z_test=data.frame(Z=Z_test)
        colnames(data_Z_test)=paste0("Z", 1:ncol(Z_test))
        
		fit=learner(O_train, data_Z_train)
        res=O_test-prediction(fit, data_Z_test)  
        Residuals[inds_test, 1]=res
	           
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
#' @param training_y Training part for Y task.
#'
.dml_association_factor_testing=function(Xrs, Y, Z, col_X, splits, factors, training_y){

  indices_x=as.matrix(1:length(col_X))
  K=length(splits)
  results <- apply(indices_x, 1, function(idx){
	  
      		i <- idx
      		
      
      		Yrs=.get_individual_residuals(Y=Y, Z=Z, ind_x=i, splits=splits, factors=factors, training_part=training_y)
      		
      		stat=sum(Xrs[,i]*Yrs, na.rm=T)
			stats_mean=mean(Xrs[,i]*Yrs, na.rm=T)
			var_S=sum((Xrs[,i]*Yrs-stats_mean)**2)
    
      		
      		p <- 2 * pnorm(-abs(stat / sqrt(var_S)), 0, 1) # compute asymptotic p-values based on normal distribution
      		
      		c(matrix1_col = col_X[i], p_value=p)
  })
  results<- as.data.frame(t(results))
  return(results)
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
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation.
#' @param training_y Training part for Y task.
#'
.dml_association_factor_testing_parallel=function(Xrs, Y, Z, col_X, splits, factors,  BPPARAM, training_y){

  indices_x=as.matrix(1:length(col_X))
  
  compute_results <- function(id){
	  
      
      		Yrs=.get_individual_residuals(Y=Y, Z=Z, ind_x=id, splits=splits, factors=factors, training_part=training_y)
      		
      		stat=sum(Xrs[,id]*Yrs, na.rm=T)
			stats_mean=mean(Xrs[,id]*Yrs, na.rm=T)
			var_S=sum((Xrs[,id]*Yrs-stats_mean)**2)
    
      		
      		p <- 2 * pnorm(-abs(stat / sqrt(var_S)), 0, 1) # compute asymptotic p-values based on normal distribution
      		
      		c(matrix1_col = col_X[id], p_value=p)
  }
  
  num_cores = BiocParallel::multicoreWorkers()
  results = BiocParallel::bplapply(
			split(seq_len(nrow(indices_x)), cut(seq_len(nrow(indices_x)), breaks=num_cores)), 
			function(ids) lapply(ids, compute_results), 
			BPPARAM = BPPARAM
  )
  results = do.call(c, results)
	
  results=do.call(rbind, results) 
  results=as.data.frame(results)
  rownames(results)=paste0(1:nrow(results))
  
  results<- as.data.frame(results)
  return(results)
}
