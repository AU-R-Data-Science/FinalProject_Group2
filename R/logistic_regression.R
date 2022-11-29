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

#' Get Confusion Matrix
#' @description This function computes Confusion Matrix for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{tp}{The \code{double} value of True Positives. Count of data points for which actual and predicted value is 1.}
#'  \item{tn}{The \code{double} value of True Negatives. Count of data points for which actual and predicted value is 0.}
#'  \item{fp}{The \code{double} value of False Positives. Count of data points for which actual value is 0 and predicted value is 1.}
#'  \item{fn}{The \code{double} value of False Negatives. Count of data points for which actual value is 1 and predicted value is 0.}
#'  \item{total}{The \code{double} value of total number of data points.}
#' }
#' @author Lydia Reedstrom
#' @export
get_confusion_matrix <- function(lr_result, cut_off = 0.5)
{
    p <- 1 / (1 + exp(-lr_result$predictors %*% lr_result$beta_estimate))
    predicted_response <- ifelse(p > cut_off, 1, 0)

    tn <- sum((lr_result$response == 0) & (predicted_response == 0))
    fp <- sum((lr_result$response == 0) & (predicted_response == 1))
    fn <- sum((lr_result$response == 1) & (predicted_response == 0))
    tp <- sum((lr_result$response == 1) & (predicted_response == 1))

    total <- tp + tn + fp + fn

    confusion_matrix <- list("tp" = tp,
                             "tn" = tn,
                             "fp" = fp,
                             "fn" = fn,
                             "total" = total)

    return(confusion_matrix)
}

#' Get Prevalence
#' @description This function computes Prevalence for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{double} value of the prevalence found.
#' @author Lydia Reedstrom
#' @export
get_prevalence <- function(lr_result, cut_off = 0.5)
{
    confusion_matrix <- get_confusion_matrix(lr_result, cut_off)
    prevalence <- (confusion_matrix$fn + confusion_matrix$tp) / confusion_matrix$total

    return(prevalence)
}

#' Get Accuracy
#' @description This function computes Accuracy for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{double} value of the accuracy found.
#' @author Lydia Reedstrom
#' @export
get_accuracy <- function(lr_result, cut_off = 0.5)
{
    confusion_matrix <- get_confusion_matrix(lr_result, cut_off)
    accuracy <- (confusion_matrix$tp + confusion_matrix$tn) / confusion_matrix$total

    return(accuracy)
}

#' Get Sensitivity
#' @description This function computes Sensitivity for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{double} value of the sensitivity found.
#' @author Lydia Reedstrom
#' @export
get_sensitivity <- function(lr_result, cut_off = 0.5)
{

    confusion_matrix <- get_confusion_matrix(lr_result, cut_off)
    sensitivity <- confusion_matrix$tp / (confusion_matrix$tp + confusion_matrix$fn)

    return(sensitivity)
}

#' Get Specificity
#' @description This function computes Specificity for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{double} value of the Specificity found.
#' @author Lydia Reedstrom
#' @export
get_specificity <- function(lr_result, cut_off = 0.5)
{
    confusion_matrix <- get_confusion_matrix(lr_result, cut_off)
    specificity <- confusion_matrix$tn / (confusion_matrix$tn + confusion_matrix$fp)

    return(specificity)
}

#' Get False Discovery Ratio
#' @description This function computes False Discovery Ratio for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{double} value of the False Discovery Ratio found.
#' @author Lydia Reedstrom
#' @export
get_false_discovery_ratio <- function(lr_result, cut_off = 0.5)
{
    confusion_matrix <- get_confusion_matrix(lr_result, cut_off)
    fdr <- confusion_matrix$fp / (confusion_matrix$fp + confusion_matrix$tp)

    return(fdr)
}

#' Get Diagnostic Odds Ratio
#' @description This function computes Diagnostic Odds Ratio for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{double} value of the Diagnostic Odds Ratio found.
#' @author Lydia Reedstrom
#' @export
get_diagnostic_odds_ratio <- function(lr_result, cut_off = 0.5)
{
    confusion_matrix <- get_confusion_matrix(lr_result, cut_off)

    fpr <- confusion_matrix$fp / (confusion_matrix$fp + confusion_matrix$tn)
    fnr <- confusion_matrix$fn / (confusion_matrix$fn + confusion_matrix$tp)

    sensitivity <- get_sensitivity(lr_result, cut_off)
    specificity <- get_specificity(lr_result, cut_off)

    dor <- (sensitivity / fpr) / (fnr / specificity)

    return(dor)
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
