% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{.get_residuals}
\alias{.get_residuals}
\title{Covariate adjustment and residual computation}
\usage{
.get_residuals(
  Outcome,
  Z,
  splits,
  training_part = 1,
  learner = lm_learner_simple,
  prediction = lm_predict_simple,
  transform_functions = NULL
)
}
\arguments{
\item{Outcome}{A matrix containing the outcome measurements.}

\item{Z}{A matrix containing the covariates Z.}

\item{splits}{The splits of the samples determining the K-fold cross fitting.}

\item{training_part}{The part of the training data that is used to train the prediction model. Default is 1.}

\item{learner}{Function pointer to prediction model (training). Default is 'lm_learner_simple'.}

\item{prediction}{Function pointer to prediction model (prediction). Default is 'lm_predict_simple'.}

\item{transform_functions}{List of functions describing the transformation that should be applied to the Outcome measurements. Default is NULL, leading to identity transformation.}
}
\description{
This function performs the covariate adjustment.
}
