% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic_regression.R
\name{logistic_regression}
\alias{logistic_regression}
\title{The main Logistic Regression function}
\usage{
logistic_regression(response, predictors)
}
\arguments{
\item{response}{A \code{double} value of the vector containing the response of interest.}

\item{predictors}{An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.}
}
\value{
A \code{list} from running the get_beta_estimate function and containing the following objects:
\describe{
 \item{initial_beta}{The initial values of coefficients}
 \item{beta_estimate}{The estimated coefficients of the logistic regression}
 \item{response}{The \code{double} vector containing the response used for the estimation}
 \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation}
}
}
\description{
This function adds intercept column and uses get_beta_estimate function to find coefficients of the logistic regression
}
\author{
Saksham Goel
}
