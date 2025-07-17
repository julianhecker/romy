

#' helper function to determine if x has less or equal than 5 levels
#'
#' @param x Value vector.
#'
.is_categorical <- function(x) {
  return(length(unique(x)) <= 5)
}


###################################################################################
### Linear regression with splines

#' Training of linear regression model
#'
#' This function performs training for a linear regression model, including interaction and squared terms.
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#' @param lambda Value for lambda to ensure that Regression matrix is invertible.
#'
#' @return trained model.
#'
#' @export
lm_learner=function(Outcome, data, lambda=NULL)
{
  ### learn using non-NA data
  ii = !is.na(Outcome)
  Outcome = Outcome[ii]
  data=as.matrix(data[ii,])
  
  ### add splines terms
  new_data=data.frame(matrix(ncol = 0, nrow = length(Outcome)))
  for(i in 1:ncol(data))
  {
	  if(!.is_categorical(data[,i])){new_data=cbind(new_data, bs(data[,i], df=3))}
	  if(.is_categorical(data[,i])){new_data=cbind(new_data, data[,i])}
  } 
  for(i in 1:ncol(data))
  {
	  for(j in 1:ncol(data))
	  {
		 if(j>i)
		 {
			new_data=cbind(new_data,(data[,i]-mean(data[,i]))*(data[,j]-mean(data[,j])))
		 }
	  }
  }
  ### add intercept
  new_data=cbind(new_data, 1)
  new_data=as.matrix(new_data)
  
  XX=t(new_data) %*% new_data
  if(is.null(lambda)){
	ed=eigen(XX)
	lambda <- 0.5*min(ed$values[ed$values > 10^-4])
  }
  ### least squares + regularization
  mm = XX + lambda*diag(ncol(new_data))
  xxx= solve(mm) %*% t(new_data)
  beta=xxx %*% as.matrix(Outcome)
  
  
  return(as.numeric(beta))
}



#' Prediction using linear regression model
#'
#' This function performs prediction using a linear regression model.
#'
#' @param obj Trained linear regression model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
lm_predict=function(obj, data)
{
  ### add interaction and squared terms
  new_data=data.frame(matrix(ncol = 0, nrow = nrow(data)))
  for(i in 1:ncol(data))
  {
	  if(!.is_categorical(data[,i])){new_data=cbind(new_data, bs(data[,i], df=3))}
	  if(.is_categorical(data[,i])){new_data=cbind(new_data, data[,i])}
  } 
  for(i in 1:ncol(data))
  {
	  for(j in 1:ncol(data))
	  {
		 if(j>i)
		 {
			new_data=cbind(new_data,(data[,i]-mean(data[,i]))*(data[,j]-mean(data[,j])))
		 }
	  }
  }
  
  ### add intercept
  new_data=cbind(new_data, 1)
  
  ### predict
  new_data=as.matrix(new_data)
  pred = new_data %*% obj
  
  return(as.numeric(pred))
}
###################################################################################
### Linear regression with quadratic and interaction terms

#' Training of linear regression model
#'
#' This function performs training for a linear regression model, including interaction and squared terms.
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#' @param lambda Value for lambda to ensure that Regression matrix is invertible.
#'
#' @return trained model.
#'
#' @export
lm_learner_q=function(Outcome, data, lambda=NULL)
{
  ### learn using non-NA data
  ii = !is.na(Outcome)
  Outcome = Outcome[ii]
  data=as.matrix(data[ii,])
  
  ### add splines terms
  new_data=data.frame(matrix(ncol = 0, nrow = length(Outcome)))
  for(i in 1:ncol(data))
  {
	  if(!.is_categorical(data[,i])){new_data=cbind(new_data, data[,i]) ;new_data=cbind(new_data, (data[,i]-mean(data[,i]))**2)}
	  if(.is_categorical(data[,i])){new_data=cbind(new_data, data[,i])}
  } 
  for(i in 1:ncol(data))
  {
	  for(j in 1:ncol(data))
	  {
		 if(j>i)
		 {
			new_data=cbind(new_data,(data[,i]-mean(data[,i]))*(data[,j]-mean(data[,j])))
		 }
	  }
  }
  ### add intercept
  new_data=cbind(new_data, 1)
  new_data=as.matrix(new_data)
  
  XX=t(new_data) %*% new_data
  if(is.null(lambda)){
	ed=eigen(XX)
	lambda <- 0.5*min(ed$values[ed$values > 10^-4])
  }
  ### least squares + regularization
  mm = XX + lambda*diag(ncol(new_data))
  xxx= solve(mm) %*% t(new_data)
  beta=xxx %*% as.matrix(Outcome)
  
  
  return(as.numeric(beta))
}



#' Prediction using linear regression model
#'
#' This function performs prediction using a linear regression model.
#'
#' @param obj Trained linear regression model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
lm_predict_q=function(obj, data)
{
  ### add interaction and squared terms
  new_data=data.frame(matrix(ncol = 0, nrow = nrow(data)))
  for(i in 1:ncol(data))
  {
	  if(!.is_categorical(data[,i])){new_data=cbind(new_data, data[,i]) ;new_data=cbind(new_data, (data[,i]-mean(data[,i]))**2)}
	  if(.is_categorical(data[,i])){new_data=cbind(new_data, data[,i])}
  } 
  for(i in 1:ncol(data))
  {
	  for(j in 1:ncol(data))
	  {
		 if(j>i)
		 {
			new_data=cbind(new_data,(data[,i]-mean(data[,i]))*(data[,j]-mean(data[,j])))
		 }
	  }
  }
  
  ### add intercept
  new_data=cbind(new_data, 1)
  
  ### predict
  new_data=as.matrix(new_data)
  pred = new_data %*% obj
  
  return(as.numeric(pred))
}

###################################################################################
### simple linear regression

#' Training of simple linear regression model
#'
#' This function performs training for a simple linear regression model.
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#' @param lambda Value for lambda to ensure that Regression matrix is invertible.
#'
#' @return trained model.
#'
#' @export
lm_learner_simple=function(Outcome, data, lambda=0.01)
{
  ### learn using non-NA data
  ii = !is.na(Outcome)
  Outcome = Outcome[ii]
  data=as.matrix(data[ii,])
  new_data=data
  
  ### add intercept
  new_data=cbind(new_data, 1)
  new_data=as.matrix(new_data)
  
  
  XX=t(new_data) %*% new_data
  if(is.null(lambda)){
	ed=eigen(XX)
	lambda <- 0.5*min(ed$values[ed$values > 10^-4])
  }
  ### least squares + regularization
  mm = XX + lambda*diag(ncol(new_data))
  xxx= solve(mm) %*% t(new_data)
  beta=xxx %*% as.matrix(Outcome)
  
  
  return(as.numeric(beta))
  
}


#' Prediction using simple linear regression model
#'
#' This function performs prediction using a linear regression model.
#'
#' @param obj Trained linear regression model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
lm_predict_simple=function(obj, data)
{
  ### add intercept
  data=cbind(data, 1)
  ### predict
  data=as.matrix(data)
  pred = data %*% obj
  ### return predicted values
  return(as.numeric(pred))
}


###################################################################################
### GAM

#' Training of GAM model
#'
#' This function performs training for a GAM model.
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#'
#' @return trained model.
#'
#' @export
mlr3_gam_learner=function(Outcome, data)
{
  if (!"regr.gam" %in% mlr3::mlr_learners$keys()) {
     mlr3extralearners::install_learners("regr.gam")
  }
  data=data.frame(data=data, Outcome=Outcome)
  predictors=colnames(data)[1:(ncol(data)-1)]
  ### learn using non-NA data 
  data=data[complete.cases(data),]
  
  if(!.is_categorical(data[,predictors[1]])){ smooth_terms=paste0("s(",predictors[1],", k=3)")}
  if(.is_categorical(data[,predictors[1]])){ smooth_terms=paste0(predictors[1])}
  
  if(length(predictors)>1){
	  for(i in 2:length(predictors))
	  {
		if(!.is_categorical(data[,predictors[i]])){ smooth_terms=paste0(smooth_terms,"+","s(",predictors[i],", k=3)")}
		if(.is_categorical(data[,predictors[i]])){ smooth_terms=paste0(smooth_terms, "+",predictors[i])}
	  }
  }
  formula <- as.formula(paste("Outcome ~", smooth_terms))

  learner <- lrn("regr.gam")

  learner$param_set$values <- list(
	formula = formula,
	method = "REML"
  )
  
  task <- TaskRegr$new("prediction", backend = data, target = "Outcome")
  
  learner$train(task)
  
  return(learner)
}

#' Prediction using GAM
#'
#' This function performs prediction using a trained GAM model.
#'
#' @param obj Trained GAM model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
mlr3_gam_predict=function(obj, data)
{
  ### setup data
  data <- data.frame(data=data)
  ### predict
  prediction <- obj$predict_newdata(data)
  ### return predicted values
  return(prediction$response)
}


###################################################################################
### LightGBM


#' Training of lightgbm model
#'
#' This function performs training for a lightgbm model.
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#'
#' @return trained model.
#'
#' @export
mlr3_lightgbm_learner=function(Outcome, data)
{
  if (!"regr.lightgbm" %in% mlr3::mlr_learners$keys()) {
     mlr3extralearners::install_learners("regr.lightgbm")
  }
  data <- data.frame(Outcome=Outcome, data=data)
  ### learn using non-NA data
  data=data[complete.cases(data),]
  learner <- mlr3::lrn("regr.lightgbm")
  task <- mlr3::TaskRegr$new("prediction", backend = data, target = "Outcome")
  
  learner$train(task)
  
  return(learner)
}


#' Prediction using lightgbm
#'
#' This function performs prediction using a trained lightgbm model.
#'
#' @param obj Trained lightgbm model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
mlr3_lightgbm_predict=function(obj, data)
{
  ### setup data
  data <- data.frame(data=data)
  ### predict
  prediction <- obj$predict_newdata(data)
  ### return predicted values
  return(prediction$response)
}





###################################################################################
### SuperLearner


#' Training of SuperLearner for gaussian link functions
#'
#' This function performs training for SuperLearner approach with gaussian link function
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#'
#' @return trained model.
#'
#' @export
superlearner_gaussian_learner=function(Outcome, data)
{
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' is required for this function. Please install it.", call. = FALSE)
  }
  learners <- c("SL.glm", "SL.rpart", "SL.randomForest")
  cn=colnames(data)
  data <- data.frame(Outcome=Outcome, data=data)
  ### learn using non-NA data
  data=data[complete.cases(data),]
  X=as.data.frame(data[,-1])
  colnames(X)=cn
  fit <- SuperLearner(Y = data$Outcome, X = X, 
                         family = gaussian(), 
                         SL.library = learners,
                         method = "method.NNLS")  # Non-negative least squares
  
  return(fit)
}


#' Prediction using SuperLearner (gaussian link)
#'
#' This function performs prediction using a trained SuperLearner model with gaussian link.
#'
#' @param obj Trained SuperLearner model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
superlearner_gaussian_predict=function(obj, data)
{
  ### setup data
  cn=colnames(data)
  data <- data.frame(data=data)
  colnames(data)=cn
  ### predict
  preds <- predict(obj, newdata = data)$pred

  ### return predicted values
  return(preds)
}




#' Training of SuperLearner for binomial link functions
#'
#' This function performs training for SuperLearner approach with binomial link function
#'
#' @param Outcome Outcome data.
#' @param data Covariate data for prediction.
#'
#' @return trained model.
#'
#' @export
superlearner_binomial_learner=function(Outcome, data)
{
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' is required for this function. Please install it.", call. = FALSE)
  }
  learners <- c("SL.glm", "SL.rpart", "SL.randomForest")
  cn=colnames(data)
  data <- data.frame(Outcome=Outcome, data=data)
  ### learn using non-NA data
  data=data[complete.cases(data),]
  X=as.data.frame(data[,-1])
  colnames(X)=cn
  fit <- SuperLearner(Y = data$Outcome, X =X, 
                         family = binomial(), 
                         SL.library = learners,
                         method = "method.NNLS")  # Non-negative least squares
  
  return(fit)
}


#' Prediction using SuperLearner (binomial link)
#'
#' This function performs prediction using a trained SuperLearner model with binomial link.
#'
#' @param obj Trained SuperLearner model.
#' @param data Covariate data for prediction.
#'
#' @return predicted Outcome values.
#'
#' @export
superlearner_binomial_predict=function(obj, data)
{
  cn=colnames(data)
  data <- data.frame(data=data)
  colnames(data)=cn
  ### predict
  preds <- predict(obj, newdata = data)$pred

  ### return predicted values
  return(preds)
}
