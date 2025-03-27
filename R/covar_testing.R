
#' @importFrom stats pnorm pchisq as.formula complete.cases sd
#' @importFrom ACAT ACAT
#' @importFrom tibble add_row
#' @importFrom magrittr %>%
NULL


#' ROMY (co)variance testing
#'
#'
#' This function tests for (co)variance effects between variables in Y and X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param index_triples A dataframe with three columns containing the indices in X and Y to be tested, respectively. Default is NULL, leading to testing all possible triples.
#' @param learner_x Prediction model function to learn prediction of X given Z. Defaul is linear regression.
#' @param prediction_x Prediction model to predict X given Z. Defaul is linear regression.
#' @param learner_y_stage1 Prediction model function to learn prediction of Y given X and Z in stage 1. Defaul is linear regression.
#' @param prediction_y_stage1 Prediction model to predict Y given X and Z in stage 1. Defaul is linear regression.
#' @param learner_y_stage2 Prediction model function to learn prediction of Y in stage 2. Defaul is linear regression.
#' @param prediction_y_stage2 Prediction model to predict Y in stage 2. Defaul is linear regression.
#' @param method Method for controlling cross-fitting. Options are 'double_cf', 'single_cf, or 'no_split'. Default is 'double_cf'.
#' @param K Value for K for the K fold cross fitting. Default is K=5.
#' @param split_ratio Ratio for splitting of training data for the two prediction tasks. Default is 0.3333/0.3333/0.3333.
#' @param parallel Logic value indicating if parallel computation should be used. Default is FALSE. If TRUE, BPPARAM needs to be initialized
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation, default is NULL.
#'
#' return Dataframe containing the p-values for all tests.
#'
#' @export
covar_testing=function(Y, X, Z, index_triples=NULL,
learner_x=lm_learner, prediction_x=lm_predict, learner_y_stage1=lm_learner, prediction_y_stage1=lm_predict, 
learner_y_stage2=lm_learner, prediction_y_stage2=lm_predict, method="double_cf", K=5, split_ratio=c(0.25, 0.25, 0.25, 0.25), parallel=FALSE, BPPARAM=NULL)
{
	### initial checks
	.input_checks(Y=Y, X=X, Z=Z, nomissX=TRUE, nomissY=FALSE, parallel=parallel, BPPARAM=BPPARAM)
	
	.check_model(learner=learner_x, prediction=prediction_x)
	.check_model(learner=learner_y_stage1, prediction=prediction_y_stage1)
	.check_model(learner=learner_y_stage2, prediction=prediction_y_stage2)
	
	if(!(method %in% c("double_cf","single_cf","no_split")))
	{
		 stop("method unknown.")
	}
	################################################
	### get dimensions
	N=nrow(Y)
	n_X=ncol(X)
	n_Y=ncol(Y)
	if(is.null(index_triples)){ 
		index_triples <- expand.grid(1:n_X, 1:n_Y, 1:n_Y)
		index_triples_sorted <- index_triples
		index_triples_sorted[, 2:3] <- t(apply(index_triples_sorted[, 2:3], 1, sort))
		index_triples <- index_triples_sorted[!duplicated(index_triples_sorted), ]
	}
	
	################################################################
	### scale data, missing data? 
	Y=scale(Y, center=TRUE, scale = FALSE)
	X=scale(X, center=TRUE, scale = TRUE)
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
		splits=.create_splits_4subs(K=K, N=N, split_ratio=split_ratio)
	}
	if(method=="single_cf")
	{
		splits=.create_splits_4subs(K=K, N=N, split_ratio=split_ratio, single=TRUE)
	}
	if(method=="no_split")
	{
		splits=.create_no_split_data(N=N, subs=4)
	}
    ################################################################
	### compute residuals
	if(!parallel){
		Xrs=.get_residuals(Outcome=X, Z=Z, splits=splits, training_part=1, learner=learner_x, prediction=prediction_x)
	}
	if(parallel){
		Xrs=.get_residuals_parallel(Outcome=X, Z=Z, splits=splits, BPPARAM=BPPARAM, training_part=1, learner=learner_x, prediction=prediction_x)
	}
	
    ################################################################
	### perform double machine learning (co)variance testing
	if(!parallel){
		results_df=.covar_testing_sub(X=X, Xrs=Xrs, Y=Y, col_X=colnames(X), col_Y=colnames(Y), Z=Z, K=K, splits=splits, index_triples=index_triples, 
		learner_x=learner_x, prediction_x=prediction_x, learner_y_stage1=learner_y_stage1, learner_y_stage2=learner_y_stage2, 
		prediction_y_stage1=prediction_y_stage1, prediction_y_stage2=prediction_y_stage2)
	}
	if(parallel){
		results_df=.covar_testing_sub_parallel(X=X, Xrs=Xrs, Y=Y, col_X=colnames(X), col_Y=colnames(Y), Z=Z, K=K, splits=splits, BPPARAM=BPPARAM, index_triples=index_triples, 
		learner_x=learner_x, prediction_x=prediction_x, learner_y_stage1=learner_y_stage1, learner_y_stage2=learner_y_stage2, 
		prediction_y_stage1=prediction_y_stage1, prediction_y_stage2=prediction_y_stage2)
	}
	
    ###################################################################
	### formatting results
    colnames(results_df)[1:3] <- c("X_id", "Y_id1", "Y_id2")
	  
    results_df[, 4] <- as.numeric(results_df[,4])
    return(results_df)
      
}


#' This function computes the test statistic for (co)variance effects
#'
#' @param X A matrix containing the measurements for X.
#' @param Xrs Matrices containing the residuals for X.
#' @param Y A matrix containing the measurements for Y. 
#' @param Z A matrix containing the measurements for Z. 
#' @param col_X Column names for X.
#' @param col_Y Column names for Y.
#' @param K Number of cross-folds.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param index_triples A dataframe with three columns containing the indices in X and Y to be tested, respectively.
#' @param learner_x Prediction model function to learn prediction of X given Z. 
#' @param prediction_x Prediction model to predict X given Z. 
#' @param learner_y_stage1 Prediction model function to learn prediction of Y given X and Z in stage 1. 
#' @param prediction_y_stage1 Prediction model to predict Y given X and Z in stage 1. 
#' @param learner_y_stage2 Prediction model function to learn prediction of Y in stage 2. 
#' @param prediction_y_stage2 Prediction model to predict Y in stage 2. 
#'
.covar_testing_sub=function(X, Xrs, Y, Z, col_X, col_Y, K, splits, index_triples, learner_x, 
prediction_x, learner_y_stage1, learner_y_stage2, prediction_y_stage1, prediction_y_stage2){

	results <- apply(index_triples, 1, function(idx){
	  
			i = idx[1]
			j1 = idx[2]
			j2 = idx[3]
			
			stat=0
			var_est=0
			
			for(k in 1:K){
				
				inds_train1=splits[[k]]$inds_train[[2]]
				inds_train2=splits[[k]]$inds_train[[3]]
				inds_train3=splits[[k]]$inds_train[[4]]
				inds_test=splits[[k]]$inds_test
				
				Y1_1=Y[inds_train1,j1]; 
				Y2_2=Y[inds_train2,j2]
				Y1_3=Y[inds_train3,j1]; 
				Y2_3=Y[inds_train3,j2]
				
				X_1=X[inds_train1,i];
				X_2=X[inds_train2,i];
				X_3=X[inds_train3,i];
				
				Z_1=Z[inds_train1,];
				Z_2=Z[inds_train2,];
				Z_3=Z[inds_train3,];
				
				
				train=.train_models_covar(Y1_1=Y1_1, Y2_2=Y2_2, Y1_3=Y1_3, Y2_3=Y2_3, X_1=X_1, X_2=X_2, X_3=X_3, Z_1=Z_1, Z_2=Z_2, Z_3=Z_3,
				learner_y_stage1=learner_y_stage1, learner_y_stage2=learner_y_stage2, 
				learner_x=learner_x, prediction_y_stage1=prediction_y_stage1, variance=j1==j2)
				
				eval=.evaluate_models_covar(Y1=Y[inds_test, j1], Y2=Y[inds_test, j2], X=X[inds_test], Z=Z[inds_test,], 
				model_data=train, prediction_y_stage1=prediction_y_stage1, prediction_y_stage2=prediction_y_stage2, prediction_x=prediction_x)
				
				stat=stat + sum(Xrs[[1]][inds_test,i]*eval$res_yy, na.rm=T)
				var_est=var_est + sum(Xrs[[1]][inds_test,i]**2*eval$res_yy**2, na.rm=T)
				
			}
			zsc=stat/sqrt(var_est)
			pval=2*pnorm(-abs(zsc),0,1)
			c(matrix1_col = col_X[i], matrix2_col = col_Y[j1], matrix3_col = col_Y[j2], p_value=pval)
	  })
	results<- as.data.frame(t(results))
	return(results)
}



#' This function computes the test statistic for (co)variance effects, using parallel computation
#'
#' @param X A matrix containing the measurements for X.
#' @param Xrs Matrices containing the residuals for X.
#' @param Y A matrix containing the measurements for Y. 
#' @param Z A matrix containing the measurements for Z. 
#' @param col_X Column names for X.
#' @param col_Y Column names for Y.
#' @param K Number of cross-folds.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation.
#' @param index_triples A dataframe with three columns containing the indices in X and Y to be tested, respectively.
#' @param learner_x Prediction model function to learn prediction of X given Z.
#' @param prediction_x Prediction model to predict X given Z.
#' @param learner_y_stage1 Prediction model function to learn prediction of Y given X and Z in stage 1. D
#' @param prediction_y_stage1 Prediction model to predict Y given X and Z in stage 1. 
#' @param learner_y_stage2 Prediction model function to learn prediction of Y in stage 2. 
#' @param prediction_y_stage2 Prediction model to predict Y in stage 2. 
#'
.covar_testing_sub_parallel=function(X, Xrs, Y, Z, col_X, col_Y, K, splits, BPPARAM, index_triples, learner_x, 
prediction_x, learner_y_stage1, learner_y_stage2, prediction_y_stage1, prediction_y_stage2){

	compute_results <- function(id){
	  
			i = index_triples[id,1]
			j1 = index_triples[id,2]
			j2 = index_triples[id,3]
			
			stat=0
			var_est=0
			
			for(k in 1:K){
				
				inds_train1=splits[[k]]$inds_train[[2]]
				inds_train2=splits[[k]]$inds_train[[3]]
				inds_train3=splits[[k]]$inds_train[[4]]
				inds_test=splits[[k]]$inds_test
				
				Y1_1=Y[inds_train1,j1]; 
				Y2_2=Y[inds_train2,j2]
				Y1_3=Y[inds_train3,j1]; Y2_3=Y[inds_train3,j2]
				
				X_1=X[inds_train1,i];
				X_2=X[inds_train2,i];
				X_3=X[inds_train3,i];
				
				Z_1=Z[inds_train1,];
				Z_2=Z[inds_train2,];
				Z_3=Z[inds_train3,];
				
				
				train=.train_models_covar(Y1_1=Y1_1, Y2_2=Y2_2, Y1_3=Y1_3, Y2_3=Y2_3, X_1=X_1, X_2=X_2, X_3=X_3, Z_1=Z_1, Z_2=Z_2, Z_3=Z_3,
				learner_y_stage1=learner_y_stage1, learner_y_stage2=learner_y_stage2, 
				learner_x=learner_x, prediction_y_stage1=prediction_y_stage1, variance=j1==j2)
				
				eval=.evaluate_models_covar(Y1=Y[inds_test, j1], Y2=Y[inds_test, j2], X=X[inds_test], Z=Z[inds_test,], 
				model_data=train, prediction_y_stage1=prediction_y_stage1, prediction_y_stage2=prediction_y_stage2, prediction_x=prediction_x)
				
				stat=stat + sum(Xrs[[1]][inds_test,i]*eval$res_yy, na.rm=T)
				var_est=var_est + sum(Xrs[[1]][inds_test,i]**2*eval$res_yy**2, na.rm=T)
				
			}
			zsc=stat/sqrt(var_est)
			pval=2*pnorm(-abs(zsc),0,1)
			c(matrix1_col = col_X[i], matrix2_col = col_Y[j1], matrix3_col = col_Y[j2], p_value=pval)
	  }
	  
	#results_list <- BiocParallel::bplapply(seq_len(nrow(index_triples)), compute_results, BPPARAM = BPPARAM)
	num_cores = BiocParallel::multicoreWorkers()
	results = BiocParallel::bplapply(
			split(seq_len(nrow(index_triples)), cut(seq_len(nrow(index_triples)), breaks=num_cores)), 
			function(ids) lapply(ids, compute_results), 
			BPPARAM = BPPARAM
	  )
	results = do.call(c, results)
	
	results=do.call(rbind, results)
	results=as.data.frame(results)
	rownames(results)=paste0(1:nrow(results))
	return(results)
}





#' This function performs the training for (co)variance testing.
#'
#' @param Y1_1 Measurements for Y1_1.
#' @param Y2_2 Measurements for Y2_2.
#' @param Y1_3 Measurements for Y1_3.
#' @param Y2_3 Measurements for Y2_3.
#' @param X_1 Measurements for X_1.
#' @param X_2 Measurements for X_2.
#' @param X_3 Measurements for X_3.
#' @param Z_1 Measurements for Z_1.
#' @param Z_2 Measurements for Z_2.
#' @param Z_3 Measurements for Z_3.
#' @param learner_y_stage1 Prediction model function to learn prediction of Y given X and Z in stage 1. 
#' @param learner_y_stage2 Prediction model function to learn prediction of Y given X and Z in stage 2.
#' @param learner_x Prediction model function to learn prediction of X given Z.
#' @param prediction_y_stage1 Prediction model to predict Y in stage 1. 
#' @param variance Logical variable to indicate if Y1=Y2. 
#'
#'
.train_models_covar=function(Y1_1, Y2_2, Y1_3, Y2_3, X_1, X_2, X_3, Z_1, Z_2, Z_3,
learner_y_stage1, learner_y_stage2, learner_x, prediction_y_stage1, variance=FALSE)
{
  
	  data_XZ1=data.frame(X=X_1, Z=Z_1)
	  data_XZ2=data.frame(X=X_2, Z=Z_2)
	  data_XZ3=data.frame(X=X_3, Z=Z_3)
	  data_Z3=data.frame(Z=Z_3)
	  
	  colnames(data_XZ1)=c(paste0("X",1),paste0("Z",1:ncol(Z_1)))
	  colnames(data_XZ2)=c(paste0("X",1),paste0("Z",1:ncol(Z_2)))
	  colnames(data_Z3)=paste0("Z",1:ncol(Z_3))
	  
	  fit1=learner_y_stage1(Y1_1, data_XZ1)
	  fit2=learner_y_stage1(Y2_2, data_XZ2)
	  
	  tmp=(Y1_3-prediction_y_stage1(fit1, data_XZ3))*(Y2_3-prediction_y_stage1(fit2, data_XZ3))
	  
	  fit=learner_y_stage2(tmp, data=data_Z3)
	  
	  
	  return(list(fit1=fit1, fit2=fit2, fit=fit))
}

#' This function evaluates the trained prediction models for (co)variance testing.
#'
#' @param Y1 Measurements for Y1.
#' @param Y2 Measurements for Y2.
#' @param X Measurements for X.
#' @param Z Measurements for Z.
#' @param model_data Fitted models.
#' @param prediction_y_stage1 Prediction model to predict Y in stage 1. 
#' @param prediction_y_stage2 Prediction model to predict Y in stage 2. 
#' @param prediction_x Prediction model to predict X. 
#' @param variance Logical variable to indicate if Y1=Y2. 
#'
.evaluate_models_covar=function(Y1, Y2, X, Z, model_data, prediction_y_stage1, prediction_y_stage2, prediction_x, variance)
{
  
	  data_XZ=data.frame(X=X, Z=Z)
	  colnames(data_XZ)=c(paste0("X",1),paste0("Z",1:ncol(Z)))
	  data_Z=data.frame(Z=Z)
	  colnames(data_Z)=paste0("Z",1:ncol(Z))
	  
	  res_y1=(Y1-prediction_y_stage1(model_data$fit1, data_XZ))
	  res_y2=(Y2-prediction_y_stage1(model_data$fit2, data_XZ))
	  
	  res_yy=res_y1*res_y2-prediction_y_stage2(model_data$fit, data_Z)
	  
	  return(list(res_y1=res_y1, res_y2=res_y2, res_yy=res_yy))
}