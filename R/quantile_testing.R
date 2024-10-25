#' ROMY association testing using quantiles
#'
#' This function performs association testing between variables in Y and X, adjusting for Z.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param Z A matrix containing the measurements for Z. No missing values allowed.
#' @param index_pairs A dataframe with two columns containing the pairs of indices in X and Y to be tested, respectively. Default is NULL, leading to testing all possible pairs.
#' @param learner_x Prediction model function to learn prediction of X given Z. Default is linear regression.
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
association_testing_quantile=function(Y, X, Z, index_pairs=NULL, 
learner_x=lm_learner, prediction_x=lm_predict, learner_y=lm_learner, prediction_y=lm_predict, 
K=5, split_ratio=c(0.5,0.5), parallel=FALSE, BPPARAM=NULL)
{
	  ### initial checks
	  .input_checks(Y=Y, X=X, Z=Z, split_ratio=split_ratio, num_splits=2, nomissX=FALSE, nomissY=TRUE, K=K, parallel=parallel, BPPARAM=BPPARAM)
	  .check_model(learner=learner_x, prediction=prediction_x)
	  
	  
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
	  fun <- function(x) { return(x) }
	  transform_functions_x=list()
	  transform_functions_x[[1]]=fun
	  
	  if(!parallel){
		  Xrs=.get_residuals(Outcome=X, Z=Z, transform_functions=transform_functions_x, splits=splits, training_part=1, learner=learner_x, prediction=prediction_x)
		  Yrs=.get_quantile_residuals(Outcome=Y, Z=Z, splits=splits, training_part=2)

	  }
	  if(parallel){
		  Xrs=.get_residuals_parallel(Outcome=X, Z=Z, transform_functions=transform_functions_x, splits=splits, BPPARAM=BPPARAM, training_part=1, learner=learner_x, prediction=prediction_x)
		  Yrs=.get_quantile_residuals_parallel(Outcome=Y, Z=Z, splits=splits, training_part=2, BPPARAM=BPPARAM)
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





#' Quantile function and residual computation using parallel computations
#'
#' @param Outcome A matrix containing the outcome measurements.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 2.
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation, default is NULL.
#'
.get_quantile_residuals_parallel=function(Outcome, Z, splits, training_part=2, BPPARAM)
{
	  ### define weights and pre-selected quantiles
	  ws=rep(0,3)
	  tau_s=rep(0, 3)
	  
	  ws[1]=sqrt(5/18)
	  ws[2]=sqrt(4/9)
	  ws[3]=sqrt(5/18)
	  
	  tau_s[1]=(-sqrt(3/5)+1)/2
	  tau_s[2]=0.5
	  tau_s[3]=(sqrt(3/5)+1)/2
	  
	  ##############################
		
	  K=length(splits)
	  
	  compute_residuals <- function(i) {
		Residuals <- vector("list", 3)
		
		for (j in 1:3) {
			Resid <- matrix(0, nrow = nrow(Outcome), ncol = 1)
		  
			  for (k in 1:K) {
						inds_train=splits[[k]]$inds_train[[training_part]] 
						inds_test=splits[[k]]$inds_test
						
						O_train=Outcome[inds_train,i];
						Z_train=Z[inds_train,];
						
						O_test=Outcome[inds_test,i];
						Z_test=Z[inds_test,];
						
						
						
						data_Z_train=data.frame(Z=Z_train)
						colnames(data_Z_train)=paste0("Z",1:ncol(Z_train))
						
						data_Z_test=data.frame(Z=Z_test)
						colnames(data_Z_test)=paste0("Z",1:ncol(Z_test))
						
						
						fit=.train_quantile_function(Y=O_train, Z=data_Z_train, tau_value=tau_s[j])
						
						res=.get_quantile_res(fit=fit, Y=O_test, Z=data_Z_test, tau_value=tau_s[j])
						
						Resid[inds_test, 1]=res
			  }
		  
			  # Store the residuals matrix for each transformation
			  Residuals[[j]] <- Resid*ws[j]
		
		}
	    list(column = i, residuals = Residuals)
	  }
	  
	  
	  num_cores = BiocParallel::multicoreWorkers()
	  Residuals_list = BiocParallel::bplapply(
			split(seq_len(ncol(Outcome)), cut(seq_len(ncol(Outcome)), breaks=num_cores)), 
			function(ids) lapply(ids, compute_residuals), 
			BPPARAM = BPPARAM
	  )
	  Residuals_list = do.call(c, Residuals_list)
	  
	  
	
	  #Residuals_list <- BiocParallel::bplapply(seq_len(ncol(Outcome)), compute_residuals, BPPARAM = BPPARAM)

	  # Reorder results based on the 'column' field to ensure correct column order
	  Residuals_list = Residuals_list[order(sapply(Residuals_list, function(x) x$column))]

	  # Combine the residuals into matrices, maintaining order
	  Residuals <- vector("list", 3)
	  for (j in 1:3) {
		
		ordered_residuals= Residuals_list |> 
		  lapply(function(x) x$residuals[[j]]) |> 
		  do.call(what = "cbind")
		
		Residuals[[j]]= ordered_residuals
	  }

	  return(Residuals)
	  
}

#' Quantile function and residual computation
#'
#' @param Outcome A matrix containing the outcome measurements.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 2.
#'
.get_quantile_residuals=function(Outcome, Z, splits, training_part=2)
{
	  ### define weights and pre-selected quantiles
	  ws=rep(0,3)
	  tau_s=rep(0, 3)
	  
	  ws[1]=sqrt(5/18)
	  ws[2]=sqrt(4/9)
	  ws[3]=sqrt(5/18)
	  
	  tau_s[1]=(-sqrt(3/5)+1)/2
	  tau_s[2]=0.5
	  tau_s[3]=(sqrt(3/5)+1)/2
	  
	  ##############################
		
	  K=length(splits)
	  Residuals=list()
	  for(j in 1:3)
	  {
			Resid=matrix(0, nrow=nrow(Outcome),ncol=ncol(Outcome))
			for(i in 1:ncol(Outcome))
			{
				  for(k in 1:K)
				  {
						inds_train=splits[[k]]$inds_train[[training_part]] 
						inds_test=splits[[k]]$inds_test
						
						O_train=Outcome[inds_train,i];
						Z_train=Z[inds_train,];
						
						O_test=Outcome[inds_test,i];
						Z_test=Z[inds_test,];
						
						
						
						data_Z_train=data.frame(Z=Z_train)
						colnames(data_Z_train)=paste0("Z",1:ncol(Z_train))
						
						data_Z_test=data.frame(Z=Z_test)
						colnames(data_Z_test)=paste0("Z",1:ncol(Z_test))
						
						
						fit=.train_quantile_function(Y=O_train, Z=data_Z_train, tau_value=tau_s[j])
						
						res=.get_quantile_res(fit=fit, Y=O_test, Z=data_Z_test, tau_value=tau_s[j])
						
						Resid[inds_test, i]=res
			
				  }
			}
			Residuals[[j]]=Resid*ws[j]
	  }
	  return(Residuals)
}

#' Function to train quantile regression using quantreg R package
#'
#' @param Y Outcome measurements.
#' @param Z Covariates.
#' @param tau_value Quantile value for quantile regression.
#'
.train_quantile_function=function(Y, Z, tau_value)
{
  ### add intercept
  ZZ=cbind(1, Z)
  
  fit=rq.fit.br(ZZ, Y, tau = tau_value)
  return(as.numeric(fit$coefficients))
  
}

#' Function to extract predicted residuals from quantile regression using quantreg R package
#'
#' @param fit Fitted object from quantile regression.
#' @param Y Outcome measurements.
#' @param Z Covariates.
#' @param tau_value Quantile value for quantile regression.
#'
.get_quantile_res=function(fit, Y, Z, tau_value)
{
  ### add intercept
  ZZ=cbind(1, Z)
  
  tmp=Y-as.matrix(ZZ) %*% fit
  res=tau_value-ifelse(tmp<0,1,0)
  return(res)
}