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
#' @param method Method for controlling cross-fitting. Options are 'double_cf' and 'single_cf'. Default is 'double_cf'.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation, default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
romy_inter=function(Y, X, Z, index_pairs=NULL, learner_ace=lm_learner_simple, prediction_ace=lm_predict_simple, 
learner_y=lm_learner, prediction_y=lm_predict, method="double_cf", K=5, parallel=FALSE, BPPARAM=NULL)
{
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, Z2=Z, nomissX=TRUE, nomissY=FALSE, parallel=parallel, BPPARAM=BPPARAM)
	  
	  .check_model(learner=learner_ace, prediction=prediction_ace)
	  .check_model(learner=learner_y, prediction=prediction_y)
	  
	  if(!(method %in% c("double_cf","single_cf")))
	  {
	  	 stop("method not implemented for association testing.")
	  }
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
		splits=.create_splits(N=N, K=K, subs=3)
		training_y=1
		training_ace_inter=2
		training_ace_y=3
	  }
	  if(method=="single_cf")
	  {
		splits=.create_splits(N=N, K=K, subs=2)
		training_y=1
		training_ace_inter=2
		training_ace_y=2
	  }
	  
	  
	  ################################################################
	  ### perform interaction testing
	  
	  
	  if(!parallel){
		  results <- apply(index_pairs, 1, function(idx){
				i <- idx[1]
				j <- idx[2]
				Y_obj=.get_residuals_and_predictions(Outcome=Y, Z=cbind(X[,i],Z), splits=splits, training_part=training_y, learner=learner_y, prediction=prediction_y)
				#Y_resid=.get_linear_additive_residuals(Outcome=Y, X=as.matrix(X[,i]), Z=Z, splits=splits,  training_part=training_y, learner=learner_y, prediction=prediction_y)
				interaction_term=X[,i]*Z[,j]
				interaction_variable_resid=.ace(Outcome=interaction_term, X=as.matrix(X[,i]), Z1=as.matrix(Z[,j]), Z2=as.matrix(Z[,-j]), splits, training_part=training_ace_inter, learner=learner_ace,
				prediction=prediction_ace)
				Y_pred=Y_obj$Predictions[[1]][,1]
				P_Y_pred=.ace(Outcome=Y_pred, X=as.matrix(X[,i]), Z1=as.matrix(Z[,j]), Z2=as.matrix(Z[,-j]), splits, training_part=training_ace_y, learner=learner_ace,
				prediction=prediction_ace)
				Y_resid=Y_obj$Residuals[[1]][,1]
				Y_resid=Y_resid+P_Y_pred
				
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
				Y_obj=.get_residuals_and_predictions(Outcome=Y, Z=cbind(X[,i],Z), splits=splits, training_part=training_y, learner=learner_y, prediction=prediction_y)
				#Y_resid=.get_linear_additive_residuals(Outcome=Y, X=as.matrix(X[,i]), Z=Z, splits=splits,  training_part=training_y, learner=learner_y, prediction=prediction_y)
				interaction_term=X[,i]*Z[,j]
				interaction_variable_resid=.ace(Outcome=interaction_term, X=as.matrix(X[,i]), Z1=as.matrix(Z[,j]), Z2=as.matrix(Z[,-j]), splits, training_part=training_ace_inter, learner=learner_ace,
				prediction=prediction_ace)
				Y_pred=Y_obj$Predictions[[1]][,1]
				P_Y_pred=.ace(Outcome=Y_pred, X=as.matrix(X[,i]), Z1=as.matrix(Z[,j]), Z2=as.matrix(Z[,-j]), splits, training_part=training_ace_y, learner=learner_ace,
				prediction=prediction_ace)
				Y_resid=Y_obj$Residuals[[1]][,1]
				Y_resid=Y_resid+P_Y_pred
				
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
#' @param Outcome Outcome measurements.
#' @param X X observations (vector).
#' @param Z1 Z1 observations (vector).
#' @param Z2 A matrix containing the covariates Z2.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Prediction model to predict Y given X and Z.
#'
.ace=function(Outcome, X, Z1, Z2, splits, training_part, learner, prediction)
{
	K=length(splits)
	resid_variable=Outcome
	for(k in 1:K){
	
			inds_train=splits[[k]]$inds_train[[training_part]] 
			inds_test=splits[[k]]$inds_test
			
			out_train=Outcome[inds_train];
			Z_train=cbind(Z1[inds_train,], Z2[inds_train,]);
			X_train=cbind(X[inds_train,], Z2[inds_train,]);
			
			out_test=Outcome[inds_test];
			Z_test=cbind(Z1[inds_test,], Z2[inds_test,]);
			X_test=cbind(X[inds_test,], Z2[inds_test,]);
			
			
			ctr=1
			while(ctr<5)
			{
			   fit=learner(out_train, X_train)
			   out_train=out_train - prediction(fit, X_train)
			   #
			   out_test=out_test - prediction(fit, X_test)
			   ##############################################################
			   fit=learner(out_train, Z_train)
			   out_train=out_train - prediction(fit, Z_train)
			   #
			   out_test=out_test - prediction(fit, Z_test)
			   ctr=ctr+1
			}
			resid_variable[inds_test]=out_test
	}
	return(resid_variable)
}

