get_initial_beta <- function(y, X)
{
    solve(t(X) %*% X) %*% t(X) %*% y
}

logistic_regression <- function(beta, y, X)
{
    p <- 1 / (1 + exp(-X %*% beta))
    sum(-((1 - y) * log(1 - p)) - (y * log(p)))
}

#' Estimate Logistic model via optimization
#' @description This function computes estimated coefficients of the logistic regression
#' @param response A \code{double} value of the vector containing the response of interest.
#' @param predictors An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression}
#'  \item{response}{The \code{double} vector containing the response used for the estimation}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation}
#' }
#' @author Parker Randall Elliott
#' @export
get_beta_estimate <- function(response, predictors)
{
    initial_beta <- get_initial_beta(response, predictors)
    beta_est <- optim(get_initial_beta(response, predictors),
                      logistic_regression,
                      y = response,
                      X = predictors)$par

    output <- list("initial_beta" = initial_beta,
                   "beta_estimate" = beta_est,
                   "response" = response,
                   "predictors" = predictors)
    class(output) <- "logistic_regression"

    return(output)
}

#' Finding Confidence Intervals
#' @description This function computes and prints bootstrap confidence intervals for each coefficient found with the 'get_beta_estimate' function.
#' @param response A \code{double} value of the vector containing the response of interest.
#' @param predictors An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param alpha A \code{double} value of the alpha used to find in the quantile function to find the confidence intervals.
#' @param replications A \code{double} value of the number of replications, or iterations, used in bootstrap to find the confidence intervals.
#' @author Parker Randall Elliott
#' @export
get_beta_confidence_intervals <- function(response, predictors, alpha, replications = 20)
{
    beta_mean <- matrix(NA, ncol = ncol(predictors), nrow = replications)

    sample_X <- cbind(response, predictors)

    for(i in 1:replications)
    {
        new_sample <- sample_X[sample(1:nrow(sample_X), replace = TRUE), ]

        beta_mean[i,] <- get_beta_estimate(response = new_sample[ ,1],
                                           predictors = new_sample[,2:ncol(new_sample)])$beta_estimate
    }

    cat("Bootstrap Confidence Intervals with alpha =",
                      alpha, "and replications =", replications, ":\n")

    for(i in 1:ncol(beta_mean))
    {
        cat("Variable:", i, "\n")
        print(quantile(beta_mean[,i], c(alpha, 1 - alpha)))
    }
}

#' Simulate Dummy Data
#' @description This function simulates dummy data that can be used to test the functionality of this package.
#' @param n A \code{double} value of the size of dataset we want to create.
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @author Parker Randall Elliott
#' @export
simulate_testing_data <- function(n)
{
    set.seed(1)
    Intercept <- rep(1, n)
    gender <- sample(c(0, 1), size = n, replace = TRUE)
    age <- round(runif(n, 18, 80))

    predictors <- cbind(Intercept, gender, age)

    xb <- -9 + 3.5*gender + 0.2*age
    p <- 1/(1 + exp(-xb))
    response <- rbinom(n = n, size = 1, prob = p)

    return(list("response" = response,
                "predictors" = predictors))
}
