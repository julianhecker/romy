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
#' @param transform_functions List of functions describing the transformation that should be applied to the Outcome measurements. Default is NULL, leading to b(x)=x.
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
						Z_train=as.matrix(Z[inds_train,]);
						
						O_test=Outcome[inds_test,i];
						Z_test=as.matrix(Z[inds_test,]);
						
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

#' Covariate adjustment and residual computation, including predictions as output
#'
#' This function performs the covariate adjustment.
#'
#' @param Outcome A matrix containing the outcome measurements.
#' @param Z A matrix containing the covariates Z.
#' @param splits The splits of the samples determining the K-fold cross fitting.
#' @param training_part The part of the training data that is used to train the prediction model. Default is 1.
#' @param learner Function pointer to prediction model (training). Default is 'lm_learner_simple'.
#' @param prediction Function pointer to prediction model (prediction). Default is 'lm_predict_simple'.
#' @param transform_functions List of functions describing the transformation that should be applied to the Outcome measurements. Default is NULL, leading to b(x)=x.
#'
.get_residuals_and_predictions=function(Outcome, Z, splits, training_part=1, learner=lm_learner_simple, prediction=lm_predict_simple, transform_functions=NULL)
{
	  if(is.null(transform_functions))
	  {
	     fun <- function(x) { return(x) }
		 transform_functions=list()
		 transform_functions[[1]]=fun
	  }
	  
	  K=length(splits)
	  Residuals=list()
	  Preds=list()
	  for(j in 1:length(transform_functions))
	  {
			Resid=matrix(0, nrow=nrow(Outcome),ncol=ncol(Outcome))
			Pred=matrix(0, nrow=nrow(Outcome),ncol=ncol(Outcome))
			for(i in 1:ncol(Outcome))
			{
				  for(k in 1:K)
				  {
						inds_train=splits[[k]]$inds_train[[training_part]] 
						inds_test=splits[[k]]$inds_test
						
						O_train=Outcome[inds_train,i];
						Z_train=as.matrix(Z[inds_train,]);
						
						O_test=Outcome[inds_test,i];
						Z_test=as.matrix(Z[inds_test,]);
						
						O_train=transform_functions[[j]](O_train)
						O_test=transform_functions[[j]](O_test)
						
						data_Z_train=data.frame(Z=Z_train)
						colnames(data_Z_train)=paste0("Z",1:ncol(Z_train))
						
						data_Z_test=data.frame(Z=Z_test)
						colnames(data_Z_test)=paste0("Z",1:ncol(Z_test))
						
						fit=learner(O_train, data_Z_train)
						
						pred=prediction(fit, data_Z_test)
						res=O_test-pred
						
						Resid[inds_test, i]=res
						Pred[inds_test, i]=pred
			
				  }
			}
			Residuals[[j]]=Resid
			Preds[[j]]=Pred
	  }
	  return(list(Residuals=Residuals, Predictions=Preds))
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
				Z_train <- as.matrix(Z[inds_train, ])
				
				O_test <- Outcome[inds_test, i]
				Z_test <- as.matrix(Z[inds_test, ])
				
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

#' K-fold 
#'
#' This function creates a partitioning of the samples for K-fold cross fitting, using subs sets for training.
#'
#' @param N Overall sample size.
#' @param K Number of folds.
#' @param subs Number of sub samples.
#'
.create_splits=function(N, K, subs)
{

    if(K==1)
    {
       indices=.create_no_split_data(N=N)
       return(indices)
    }
  
    indices=list()
    inds=1:N
    
    
    folds=sample(cut(seq(1,N), breaks=K, labels=FALSE))
    for(k in 1:K)
    {
      tmp=inds[folds==k] # test data
      tmpc=inds[folds!=k] # training data
      splits=sample(1:subs, size=length(tmpc), replace=TRUE) # split training data
      inds_train=list()
      for(s in 1:subs)
      {
        inds_train[[s]]=tmpc[splits==s]; 
      }
      tmp=list(inds_train=inds_train, inds_test=tmp)
      indices[[k]]=tmp
    }
    return(indices)
}







#' no-cross-fitting-approach
#'
#' This function sets up the data such that evaluation and training data is the same.
#'
#' @param N Overall sample size.
#'
.create_no_split_data=function(N)
{
	indices=list()
	inds=1:N
	inds_train=list()
	
	inds_train[[1]]=inds
	
	tmp=list(inds_train=inds_train,inds_test=inds)
	indices[[1]]=tmp
	return(indices)
}





#' Input check
#'
#' This function checks the input provided.
#'
#' @param Y Matrix Y.
#' @param X Matrix X.
#' @param Z Matrix Z.
#' @param Z2 Matrix Z2.
#' @param nomissX Logic value if missing values in X should be allowed.
#' @param nomissY Logic value if missing values in Y should be allowed.
#' @param parallel Logic value to indicate if parallel computing should be used.
#' @param BPPARAM BPPARAM object for BiocParallel parallel computation.
#'
.input_checks=function(Y, X, Z, Z2, nomissX, nomissY, parallel, BPPARAM){

	###############################################
	## Y, X, and Z need to be matrices
	if(!is.matrix(Y)){stop("Y is not a matrix.")}
	if(!is.matrix(X)){stop("X is not a matrix.")}
	if(!is.matrix(Z)){stop("Z is not a matrix.")}
	if(!is.matrix(Z2)){stop("Z2 is not NULL and not a matrix.")}

	################################################
	## Missingness in data
	if(sum(is.na(Z))>0){
		stop("No missing values in Z allowed.")
	}
	if(sum(is.na(Z2))>0){
		stop("No missing values in Z2 allowed if provided.")
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
	if(!all(apply(Z2, 2, sd, na.rm=T)>0)){stop("Z2 contains constant columns.")}
	
	
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
		X=data.frame(X=X)
		### fit prediction model
		fit=learner(Outcome=Y, data=X)
		### simulate new data
		X=as.matrix(rnorm(N))
		Y=X+rnorm(N)*0.1
		X=data.frame(X=X)
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


#' Matrix product with potential NA entries
#'
#' This function computes a matrix product that can deal with NA data.
#'
#' @param A Matrix A.
#' @param B Matrix B.
#'
.na_safe_matrix_product= function(A, B) {
 
  n <- nrow(A)
  m <- ncol(B)
  
  
  result <- matrix(NA_real_, nrow = n, ncol = m)
  
  for (i in 1:n) {
    for (j in 1:m) {
      a_row <- A[i, ]
      b_col <- B[, j]
      
      result[i, j] <- sum(a_row * b_col, na.rm=T)
    }
  }
  
  return(result)
}