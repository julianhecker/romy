#' ROMY interaction testing
#'
#'
#' This function performs interaction testing.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param index_pairs A dataframe with two columns containing the pairs of indices in X and Z to be tested for interaction, respectively. Default is NULL, leading to testing all possible pairs.
#' @param learner_ace Prediction model function to learn prediction in the ACE algorithm. Defaul is linear regression.
#' @param prediction_ace Prediction model to predict in the ACE algorithm. Defaul is linear regression.
#' @param learner_y Prediction model function to learn prediction of Y given X and Z (without interactions). Defaul is linear regression.
#' @param prediction_y Prediction model to predict Y given X and Z. Defaul is linear regression.
#' @param method Method for controlling cross-fitting. Options are 'double_cf', 'single_cf', or 'no_split'. Default is 'double_cf'.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param split_ratio Ratio for splitting of training data for the two prediction tasks. Default is 0.5:0.5.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation, default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
interaction_testing=function(Y, X, Z, index_pairs=NULL, learner_ace=lm_learner_simple, prediction_ace=lm_predict_simple, 
learner_y=lm_learner, prediction_y=lm_predict, method="double_cf", K=5, split_ratio=c(0.5, 0.5), parallel=FALSE, BPPARAM=NULL)
{
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, nomissX=TRUE, nomissY=FALSE, parallel=parallel, BPPARAM=BPPARAM)
	  
	  .check_model(learner=learner_ace, prediction=prediction_ace)
	  .check_model(learner=learner_y, prediction=prediction_y)
	  
	  ################################################
	  ### get dimensions
	  N=nrow(Y)
	  if(ncol(Y)>1){
		Y=as.matrix(Y[,1])
		warning("warning: Y contains more than 1 column. Extracting first column for analysis.")
	  }
	  
	  n_X=ncol(X)
	  n_Z=ncol(Z)
	  
	  if(is.null(index_pairs)){ index_pairs <- expand.grid(1:n_X, 1:n_Z)}
	  
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
	  
	  if(method=="double_cf")
	  {
		splits=.create_splits_2subs(K=K, N=N, split_ratio=split_ratio)
	  }
	  if(method=="single_cf")
	  {
		splits=.create_splits_2subs(K=K, N=N, split_ratio=split_ratio, single=TRUE)
	  }
	  if(method=="no_split")
	  {
		splits=.create_no_split_data(N=N, subs=2)
	  }
	  ################################################################
	  ### perform interaction testing
	  
	  
	  if(!parallel){
		  results <- apply(index_pairs, 1, function(idx){
				i <- idx[1]
				j <- idx[2]
				Y_resid=.get_linear_additive_residuals(Outcome=Y, X=as.matrix(X[,i]), Z=Z, splits=splits,  training_part=1, learner=learner_y, prediction=prediction_y)
				interaction_term=X[,i]*Z[,j]
				interaction_variable_resid=.ace(interaction_term=interaction_term, X=as.matrix(X[,i]), Z=Z, splits, training_part=2, learner=learner_ace,
				prediction=prediction_ace)
				
				stat=sum(Y_resid*interaction_variable_resid, na.rm=TRUE)
				var_est=sum(Y_resid**2*interaction_variable_resid**2, na.rm=TRUE)
				zsc=stat/sqrt(var_est)
				pval=2*pnorm(-abs(zsc),0,1)
				
				c(matrix1_col = colnames(X)[i], matrix2_col = colnames(Z)[j], p_value=pval)
		  })
		  results_df <- as.data.frame(t(results))
	  }
	  if(parallel){
		  compute_results=function(id){
				i <- index_pairs[id,1]
				j <- index_pairs[id,2]
				Y_resid=.get_linear_additive_residuals(Outcome=Y, X=as.matrix(X[,i]), Z=Z, splits=splits,  training_part=1, learner=learner_y, prediction=prediction_y)
				interaction_term=X[,i]*Z[,j]
				interaction_variable_resid=.ace(interaction_term=interaction_term, X=as.matrix(X[,i]), Z=Z, splits, training_part=2, learner=learner_ace,
				prediction=prediction_ace)
				
				stat=sum(Y_resid*interaction_variable_resid, na.rm=TRUE)
				var_est=sum(Y_resid**2*interaction_variable_resid**2, na.rm=TRUE)
				zsc=stat/sqrt(var_est)
				pval=2*pnorm(-abs(zsc),0,1)
				
				c(matrix1_col = colnames(X)[i], matrix2_col = colnames(Z)[j], p_value=pval)
		  }
		  num_cores = BiocParallel::multicoreWorkers()
		  results = BiocParallel::bplapply(
				split(seq_len(nrow(index_pairs)), cut(seq_len(nrow(index_pairs)), breaks=num_cores)), 
				function(ids) lapply(ids, compute_results), 
				BPPARAM = BPPARAM
		  )
		  results = do.call(c, results)
		  
		  results=do.call(rbind, results)
		  results_df=as.data.frame(results)
		  rownames(results_df)=paste0(1:nrow(results_df))
		   
	  }
	  ###################################################################
	  ### formatting results
	  
	  colnames(results_df)[1:2] <- c("X_id", "Y_id")
	  results_df[, 3] <- as.numeric(results_df[,3])
	  return(results_df)
}





#' This functions performs the alternating condition expectation algorithm, the outcome is the interaction term.
#'
#' @param interaction_term Measurements for the interaction term.
#' @param X A matrix containing X.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Prediction model to predict Y given X and Z.
#'
.ace=function(interaction_term, X, Z, splits, training_part, learner, prediction)
{
	K=length(splits)
	resid_interaction_variable=interaction_term
	for(k in 1:K){
	
			inds_train=splits[[k]]$inds_train[[training_part]] 
			inds_test=splits[[k]]$inds_test
			
			inter_train=interaction_term[inds_train];
			Z_train=Z[inds_train,];
			X_train=as.matrix(X[inds_train,]);
			
			inter_test=interaction_term[inds_test];
			Z_test=Z[inds_test,];
			X_test=as.matrix(X[inds_test,]);
			
			
			ctr=1
			while(ctr<5)
			{
			   fit=learner(inter_train, X_train)
			   inter_train=inter_train - prediction(fit, X_train)
			   #
			   inter_test=inter_test - prediction(fit, X_test)
			   ##############################################################
			   fit=learner(inter_train, Z_train)
			   inter_train=inter_train - prediction(fit, Z_train)
			   #
			   inter_test=inter_test - prediction(fit, Z_test)
			   ctr=ctr+1
			}
			resid_interaction_variable[inds_test]=inter_test
	}
	return(resid_interaction_variable)
}

#' This function computes residuals using a prediction model that does not incorporate interactions between X and Z.
#'
#' @param Outcome A matrix containing the outcome measurements.
#' @param X A matrix containing X.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Prediction model to predict Y given X and Z.
#'
.get_linear_additive_residuals=function(Outcome, X, Z, splits, training_part, learner, prediction)
{
	  K=length(splits)
	  Resid=numeric(length(Outcome))
	  for(k in 1:K)
	  {
			inds_train=splits[[k]]$inds_train[[training_part]] 
			inds_test=splits[[k]]$inds_test
			
			O_train=Outcome[inds_train];
			Z_train=Z[inds_train,];
			X_train=X[inds_train,];
			
			O_test=Outcome[inds_test];
			Z_test=Z[inds_test,];
			X_test=X[inds_test,];
			
			
			data_train=data.frame(X=X_train, Z=Z_train)
			colnames(data_train)=paste0("C",1:ncol(data_train))
			
			data_test=data.frame(X=X_test, Z=Z_test)
			colnames(data_test)=paste0("C",1:ncol(data_test))
			
			
			fit=learner(O_train, data_train)
			res=O_test-prediction(fit, data_test)
			
			Resid[inds_test]=res

	  }
	  return(Resid)
}