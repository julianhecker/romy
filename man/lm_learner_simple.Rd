% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/learner.R
\name{lm_learner_simple}
\alias{lm_learner_simple}
\title{Training of simple linear regression model}
\usage{
lm_learner_simple(Outcome, data, lambda = 0.01)
}
\arguments{
\item{Outcome}{Outcome data.}

\item{data}{Covariate data for prediction.}

\item{lambda}{Value for lambda to ensure that Regression matrix is invertible.}
}
\value{
trained model.
}
\description{
This function performs training for a simple linear regression model.
}
