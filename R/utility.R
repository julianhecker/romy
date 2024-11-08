#' Covariate adjustment and residual computation
#'
#' This function performs the covariate adjustment.
#'
#' @param Outcome A matrix containing the outcome measurements.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Function pointer to prediction model (prediction). Default is 'lm_predict_simple'.
#' @param transform_functions List of functions describing the transformation that should be applied to the Outcome measurements. Default is NULL, leading to identity transformation.
#'
.get_residuals=function(Outcome, Z, splits, training_part=1, learner=lm_learner_simple, prediction=lm_predict_simple, transform_functions=NULL)
{
	  if(is.null(transform_functions))
	  {
	     fun <- function(x) { return(x) }
		 transform_functions=list()
		 transform_functions[[1]]=fun
	  }
	  
	  K=length(splits)
	  Residuals=list()
	  for(j in 1:length(transform_functions))
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
						
						O_train=transform_functions[[j]](O_train)
						O_test=transform_functions[[j]](O_test)
						
						data_Z_train=data.frame(Z=Z_train)
						colnames(data_Z_train)=paste0("Z",1:ncol(Z_train))
						
						data_Z_test=data.frame(Z=Z_test)
						colnames(data_Z_test)=paste0("Z",1:ncol(Z_test))
						
						fit=learner(O_train, data_Z_train)
						
						res=O_test-prediction(fit, data_Z_test)
						
						Resid[inds_test, i]=res
			
				  }
			}
			Residuals[[j]]=Resid
	  }
	  return(Residuals)
}

#' Covariate adjustment and residual computation using parallel computations
#'
#' This function performs the covariate adjustment using parallel computations.
#'
#' @param Outcome A matrix containing the outcome measurements.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Function pointer to prediction model (prediction). Default is 'lm_predict_simple'.
#' @param transform_functions List of functions describing the transformation that should be applied to the Outcome measurements. Default is NULL, leading to identity transformation.
#'
.get_residuals_parallel <- function(Outcome, Z, splits, BPPARAM, training_part = 1, learner = lm_learner_simple, 
                                    prediction = lm_predict_simple, transform_functions = NULL) {
	  if (is.null(transform_functions)) {
		fun <- function(x) { return(x) }
		transform_functions <- list(fun)
	  }
	  
	  K <- length(splits)
	  
	  # Define the function for parallel computation
	  compute_residuals <- function(i) {
		Residuals <- vector("list", length(transform_functions))
		
		for (j in seq_along(transform_functions)) {
			  Resid <- matrix(0, nrow = nrow(Outcome), ncol = 1)
		  
			  for (k in 1:K) {
				inds_train <- splits[[k]]$inds_train[[training_part]]
				inds_test <- splits[[k]]$inds_test
				
				O_train <- Outcome[inds_train, i]
				Z_train <- Z[inds_train, ]
				
				O_test <- Outcome[inds_test, i]
				Z_test <- Z[inds_test, ]
				
				# Apply the transformation to the training and test data
				O_train <- transform_functions[[j]](O_train)
				O_test <- transform_functions[[j]](O_test)
				
				data_Z_train <- data.frame(Z = Z_train)
				colnames(data_Z_train) <- paste0("Z", 1:ncol(Z_train))
				
				data_Z_test <- data.frame(Z = Z_test)
				colnames(data_Z_test) <- paste0("Z", 1:ncol(Z_test))
				
				# Fit the model and compute residuals
				fit <- learner(O_train, data_Z_train)
				res <- O_test - prediction(fit, data_Z_test)
				
				# Store residuals in the corresponding locations for this column
				Resid[inds_test, 1] <- res
			  }
		  
			  # Store the residuals matrix for each transformation
			  Residuals[[j]] <- Resid
		}
		
		# Return the residuals along with the column index
		list(column = i, residuals = Residuals)
	  }
	  
	  num_cores = BiocParallel::multicoreWorkers()
	  Residuals_list = BiocParallel::bplapply(
			split(seq_len(ncol(Outcome)), cut(seq_len(ncol(Outcome)), breaks=num_cores)), 
			function(ids) lapply(ids, compute_residuals), 
			BPPARAM = BPPARAM
	  )
	  Residuals_list = do.call(c, Residuals_list)
	  
	  

	  Residuals_list = Residuals_list[order(sapply(Residuals_list, function(x) x$column))]

	  # Combine the residuals into matrices, maintaining order
	  Residuals <- vector("list", length(transform_functions))
	  for (j in seq_along(transform_functions)) {
		
		ordered_residuals= Residuals_list |> 
		  lapply(function(x) x$residuals[[j]]) |> 
		  do.call(what = "cbind")
		
		Residuals[[j]]= ordered_residuals
	  }

	  return(Residuals)
}







#' K-fold with 1 training part
#'
#' This function creates a partitioning of the samples for K-fold cross fitting, using 1 set for training.
#'
#' @param K Number of folds.
#' @param N Overall sample size.
#'
.create_splits_1sub=function(K, N){

	  indices=list()
	  inds=1:N
	  if(K>=2)
	  { 
		  folds=sample(cut(seq(1,N), breaks=K, labels=FALSE))
		  for(k in 1:K)
		  {
			tmp=inds[folds==k] # test data
			tmpc=inds[folds!=k] # training data
			
			inds_train=list()
			inds_train[[1]]=tmpc
			tmp=list(inds_train=inds_train, inds_test=tmp)
			indices[[k]]=tmp
		  }
	  }
	  if(K==1) # test data = training data = whole sample
	  {
	      inds_train=list()
		  inds_train[[1]]=inds
		  tmp=list(inds_train=inds_train, inds_test=inds)
		  indices[[1]]=tmp
	  }
	  return(indices)
}

#' K-fold with 2 training splits
#'
#' This function creates a partitioning of the samples for K-fold cross fitting, using 2 sets for training.
#'
#' @param K Number of folds.
#' @param N Overall sample size.
#' @param split_ratio Fractions describing the allocation of training data for different training tasks.
#'
.create_splits_2subs=function(K, N, split_ratio){

	  indices=list()
	  inds=1:N
	  if(K>=2)
	  { 
		  folds=sample(cut(seq(1,N), breaks=K, labels=FALSE))
		  for(k in 1:K)
		  {
			tmp=inds[folds==k] # test data
			tmpc=inds[folds!=k] # training data
			splits=sample(1:2, size=length(tmpc), replace=TRUE, prob=split_ratio) # split training data
			inds_train=list()
			inds_train[[1]]=tmpc[splits==1];inds_train[[2]]=tmpc[splits==2]
			tmp=list(inds_train=inds_train, inds_test=tmp)
			indices[[k]]=tmp
		  }
	  }
	  if(K==1) # test data = training data = whole sample
	  {
	      inds_train=list()
		  inds_train[[1]]=inds; inds_train[[2]]=inds
		  tmp=list(inds_train=inds_train, inds_test=inds)
		  indices[[1]]=tmp
	  }
	  return(indices)
}

#' K-fold with 3 training splits
#'
#' This function creates a partitioning of the samples for K-fold cross fitting, using 3 sets for training.
#'
#' @param K Number of folds.
#' @param N Overall sample size.
#' @param split_ratio Fractions describing the allocation of training data for different training tasks.
#'
.create_splits_3subs=function(K, N, split_ratio)
{
	  inds=1:N
	  indices=list()
	  if(K>=2)
	  {
		  folds=sample(cut(seq(1,N), breaks=K, labels=FALSE))
		 
		  for(k in 1:K)
		  {
			tmp=inds[folds==k] # test data
			tmpc=inds[folds!=k] # training data
			splits=sample(1:3, size=length(tmpc), replace=TRUE, prob=split_ratio)
			inds_train=list()
			inds_train[[1]]=tmpc[splits==1];inds_train[[2]]=tmpc[splits==2];inds_train[[3]]=tmpc[splits==3]
			tmp=list(inds_train=inds_train, inds_test=tmp) # split training data
			indices[[k]]=tmp
		  }
	  }
	  if(K==1) # test data = training data = whole sample
	  {
		  inds_train=list()
		  inds_train[[1]]=inds;inds_train[[2]]=inds;inds_train[[3]]=inds
		  tmp=list(inds_train=inds_train,inds_test=inds)
		  indices[[1]]=tmp
	  }
	  return(indices)
}





#' Input check
#'
#' This function checks the input provided.
#'
#' @param Y Matrix Y.
#' @param X Matrix X.
#' @param Z Matrix Z.
#' @param split_ratio Split ratio for training.
#' @param num_splits Length of split ratio vector.
#' @param nomissX Logic value if missing values in X should be allowed.
#' @param nomissY Logic value if missing values in Y should be allowed.
#' @param K K value for cross-fitting.
#' @param parallel Logic value to indicate if parallel computing should be used.
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation.
#'
.input_checks=function(Y, X, Z, split_ratio, num_splits, nomissX, nomissY, K, parallel, BPPARAM){

	###############################################
	## Y, X, and Z need to be matrices
	if(!is.matrix(Y)){stop("Y is not a matrix.")}
	if(!is.matrix(X)){stop("X is not a matrix.")}
	if(!is.matrix(Z)){stop("Z is not a matrix.")}

	################################################
	## Missingness in data
	if(sum(is.na(Z))>0){
		stop("No missing values in Z allowed.")
	}
	if(sum(is.na(X))>0 & nomissX){
		stop("No missing values in X allowed.")
	}
	if(sum(is.na(Y))>0 & nomissY){
		stop("No missing values in Y allowed.")
	}
	################################################
	if(!nomissY){if(any(apply(Y, 2, function(x){sum(is.na(x))>0.05*length(x)}))){stop("Some columns of Y contain excessive missing rates (>5%).")}}
	if(!nomissX){if(any(apply(X, 2, function(x){sum(is.na(x))>0.05*length(x)}))){stop("Some columns of X contain excessive missing rates (>5%).")}}
	
	## all measurements have to be non-constant
	if(!all(apply(Y, 2, sd, na.rm=T)>0)){stop("Y contains constant columns.")}
	if(!all(apply(X, 2, sd, na.rm=T)>0)){stop("X contains constant columns.")}
	if(!all(apply(Z, 2, sd, na.rm=T)>0)){stop("Z contains constant columns.")}
	
	
	################################################
	## Split ratio and K have to be valid
	if(!is.null(split_ratio)){
		if(any(split_ratio<0)){stop("split_ratio contains negative values.")}
		if(!(sum(split_ratio)==1)){
			stop("split_ratio doest not sum to 1.")
		}
		##
		if(!(K==1 | (K>1 & length(split_ratio)==num_splits))){stop("K or split_ratio not correctly specified.")}
	}
	################################################
	## Parallel computing objects
	if(parallel & is.null(BPPARAM)){
	   stop("Parallel computing requested but BPPARAM is not initialized.")
	}

}

#' Prediction model checker
#'
#' This function checks if the provided functions for training and prediction perform well on a simple task.
#'
#' @param learner Function pointer for prediction model (training).
#' @param prediction Function pointer for prediction model (prediction).
#'
.check_model=function(learner, prediction){
	tryCatch({
		### simulate simple model
		N=1000
		X=as.matrix(rnorm(N))
		Y=X+rnorm(N)*0.1
		### fit prediction model
		fit=learner(Outcome=Y, data=X)
		### simulate new data
		X=as.matrix(rnorm(N))
		Y=X+rnorm(N)*0.1
		### predict using trained model
		pred=prediction(fit, data=X)
		### check prediction quality, should have R^2>50% for this simple model.
		fit=lm(Y~pred)
		if(summary(fit)$r.squared<=0.5){stop("learning/prediction model does not perform well enough.")}
		return(summary(fit)$r.squared)
	}, error=function(e){
		stop("An error occured while testing the learning/prediction models: ", e$message)
	})
}



