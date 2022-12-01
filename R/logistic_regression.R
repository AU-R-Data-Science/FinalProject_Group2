get_initial_beta <- function(y, X)
{
    solve(t(X) %*% X) %*% t(X) %*% y
}

optimizing_function <- function(beta, y, X)
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
                      optimizing_function,
                      y = response,
                      X = predictors)$par

    output <- list("initial_beta" = initial_beta,
                   "beta_estimate" = beta_est)
    return(output)
}

#' The main Logistic Regression function
#' @description This function adds intercept column and uses get_beta_estimate function to find coefficients of the logistic regression
#' @param response A \code{double} value of the vector containing the response of interest.
#' @param predictors An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @return A \code{list} from running the get_beta_estimate function and containing the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression}
#'  \item{response}{The \code{double} vector containing the response used for the estimation}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation}
#' }
#' @author Saksham Goel
#' @export
logistic_regression <- function(response, predictors)
{
    Intercept <- rep(1, length(response))
    predictors <- cbind(Intercept, predictors)

    result <- get_beta_estimate(response, predictors)

    lr_result <- list("initial_beta" = result$initial_beta,
                      "beta_estimate" = result$beta_estimate,
                      "response" = response,
                      "predictors" = predictors)

    class(lr_result) <- "logistic_regression"

    return(lr_result)
}

#' Finding Confidence Intervals
#' @description This function computes and prints bootstrap confidence intervals for each coefficient found with the 'get_beta_estimate' function.
#' @param response A \code{double} value of the vector containing the response of interest.
#' @param predictors An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param alpha A \code{double} value of the alpha used to find in the quantile function to find the confidence intervals.
#' @param replications A \code{double} value of the number of replications, or iterations, used in bootstrap to find the confidence intervals.
#' @author Parker Randall Elliott
#' @export
get_beta_confidence_intervals <- function(lr_result, alpha = 0.05, replications = 20)
{
    beta_mean <- matrix(NA, ncol = ncol(lr_result$predictors), nrow = replications)

    sample_X <- cbind(lr_result$response, lr_result$predictors)

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

#' Plot Confusion Matrix
#' @description This function plots the Confusion Matrix for the found coefficients of logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @param cut_off A \code{double} value of the cut_off used to set the predicted values to 0 or 1.
#' @return A \code{list} value of the ggplot object created for the specific confusion matrix.
#' @author Saksham Goel
#' @export
plot_confusion_matrix <- function(lr_result, cut_off = 0.5)
{
    cm <- get_confusion_matrix(lr_result, cut_off)

    Target <- factor(c(0, 0, 1, 1))
    Prediction <- factor(c(0, 1, 0, 1))
    values <- c(cm$tn, cm$fp, cm$fn, cm$tp)
    df <- data.frame(Target, Prediction, values)

    plot <- ggplot2::ggplot(data =  df, mapping = ggplot2::aes(x = Target, y = Prediction)) +
        ggplot2::geom_tile(ggplot2::aes(fill = values), colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%1.0f", values)), vjust = 1) +
        ggplot2::scale_fill_gradient(low = "azure1", high = "lightblue4") +
        ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")

    print(plot)
    return(plot)
}

#' Plot Prevalence
#' @description This function plots the Prevalence for the found coefficients of logistic regression over different cut off values from 0.1 to 0.9.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{df}{The data frame containing cut-off values and the corresponding prevalence values.}
#'  \item{plot}{The ggplot object created for the prevalence bar graph.}
#' }
#' @author Saksham Goel
#' @export
plot_prevalence <- function(lr_result)
{
    prevalence_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(prevalence_df) <- c("Cut_Off", "Prevalence")

    for(i in seq(from = 0.1, to = 0.9, by = 0.1))
    {
        prevalence <- get_prevalence(lr_result, cut_off = i)
        prevalence_df[nrow(prevalence_df) + 1,] <- c(i, prevalence)
    }

    plot <- ggplot2::ggplot(data = prevalence_df, ggplot2::aes(x = Cut_Off, y = Prevalence, fill = Cut_Off)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_x_continuous("Cut Off Values",
                           labels = as.character(prevalence_df$Cut_Off),
                           breaks = prevalence_df$Cut_Off)

    print(plot)
    return(list("df" = prevalence_df,
                "plot" = plot))
}

#' Plot Accuracy
#' @description This function plots the Accuracy for the found coefficients of logistic regression over different cut off values from 0.1 to 0.9.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{df}{The data frame containing cut-off values and the corresponding accuracy values.}
#'  \item{plot}{The ggplot object created for the accuracy bar graph.}
#' }
#' @author Saksham Goel
#' @export
plot_accuracy <- function(lr_result)
{
    accuracy_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(accuracy_df) <- c("Cut_Off", "Accuracy")

    for(i in seq(from = 0.1, to = 0.9, by = 0.1))
    {
        accuracy <- get_accuracy(lr_result, cut_off = i)
        accuracy_df[nrow(accuracy_df) + 1,] <- c(i, accuracy)
    }

    plot <- ggplot2::ggplot(data = accuracy_df, ggplot2::aes(x = Cut_Off, y = Accuracy, fill = Cut_Off)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_x_continuous("Cut Off Values",
                           labels = as.character(accuracy_df$Cut_Off),
                           breaks = accuracy_df$Cut_Off)

    print(plot)
    return(list("df" = accuracy_df,
                "plot" = plot))
}

#' Plot Sensitivity
#' @description This function plots the Sensitivity for the found coefficients of logistic regression over different cut off values from 0.1 to 0.9.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{df}{The data frame containing cut-off values and the corresponding sensitivity values.}
#'  \item{plot}{The ggplot object created for the sensitivity bar graph.}
#' }
#' @author Saksham Goel
#' @export
plot_sensitivity <- function(lr_result)
{
    sensitivity_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(sensitivity_df) <- c("Cut_Off", "Sensitivity")

    for(i in seq(from = 0.1, to = 0.9, by = 0.1))
    {
        sensitivity <- get_sensitivity(lr_result, cut_off = i)
        sensitivity_df[nrow(sensitivity_df) + 1,] <- c(i, sensitivity)
    }
    plot <- ggplot2::ggplot(data = sensitivity_df, ggplot2::aes(x = Cut_Off, y = Sensitivity, fill = Cut_Off)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_x_continuous("Cut Off Values",
                           labels = as.character(sensitivity_df$Cut_Off),
                           breaks = sensitivity_df$Cut_Off)

    print(plot)
    return(list("df" = sensitivity_df,
                "plot" = plot))
}

#' Plot specificity
#' @description This function plots the specificity for the found coefficients of logistic regression over different cut off values from 0.1 to 0.9.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{df}{The data frame containing cut-off values and the corresponding specificity values.}
#'  \item{plot}{The ggplot object created for the specificity bar graph.}
#' }
#' @author Saksham Goel
#' @export
plot_specificity <- function(lr_result)
{
    specificity_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(specificity_df) <- c("Cut_Off", "Specificity")

    for(i in seq(from = 0.1, to = 0.9, by = 0.1))
    {
        specificity <- get_specificity(lr_result, cut_off = i)
        specificity_df[nrow(specificity_df) + 1,] <- c(i, specificity)
    }
    plot <- ggplot2::ggplot(data = specificity_df, ggplot2::aes(x = Cut_Off, y = Specificity, fill = Cut_Off)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_x_continuous("Cut Off Values",
                           labels = as.character(specificity_df$Cut_Off),
                           breaks = specificity_df$Cut_Off)

    print(plot)
    return(list("df" = specificity_df,
                "plot" = plot))
}

#' Plot False Discovery Ratio
#' @description This function plots the False Discovery Ratio values for the found coefficients of logistic regression over different cut off values from 0.1 to 0.9.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{df}{The data frame containing cut-off values and the corresponding False Discovery Ratio values.}
#'  \item{plot}{The ggplot object created for the False Discovery Ratio bar graph.}
#' }
#' @author Saksham Goel
#' @export
plot_false_discovery_ratio <- function(lr_result)
{
    fdr_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(fdr_df) <- c("Cut_Off", "False_Discovery_Ratio")

    for(i in seq(from = 0.1, to = 0.9, by = 0.1))
    {
        fdr <- get_false_discovery_ratio(lr_result, cut_off = i)
        fdr_df[nrow(fdr_df) + 1,] <- c(i, fdr)
    }
    plot <- ggplot2::ggplot(data = fdr_df, ggplot2::aes(x = Cut_Off, y = False_Discovery_Ratio, fill = Cut_Off)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_x_continuous("Cut Off Values",
                           labels = as.character(fdr_df$Cut_Off),
                           breaks = fdr_df$Cut_Off)

    print(plot)
    return(list("df" = fdr_df,
                "plot" = plot))
}

#' Plot Diagnostic Odds Ratio
#' @description This function plots the Diagnostic Odds Ratio for the found coefficients of logistic regression over different cut off values from 0.1 to 0.9.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @return A \code{list} value containing the following objects:
#' \describe{
#'  \item{df}{The data frame containing cut-off values and the corresponding diagnostic odds ratio values.}
#'  \item{plot}{The ggplot object created for the diagnostic odds ratio bar graph.}
#' }
#' @author Saksham Goel
#' @export
plot_diagnostic_odds_ratio <- function(lr_result)
{
    dor_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(dor_df) <- c("Cut_Off", "Diagnostic_Odds_Ratio")

    for(i in seq(from = 0.1, to = 0.9, by = 0.1))
    {
        dor <- get_diagnostic_odds_ratio(lr_result, cut_off = i)
        dor_df[nrow(dor_df) + 1,] <- c(i, dor)
    }
    plot <- ggplot2::ggplot(data = dor_df, ggplot2::aes(x = Cut_Off, y = Diagnostic_Odds_Ratio, fill = Cut_Off)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_x_continuous("Cut Off Values",
                           labels = as.character(dor_df$Cut_Off),
                           breaks = dor_df$Cut_Off)

    print(plot)
    return(list("df" = dor_df,
                "plot" = plot))
}

#' Plot Logistic curve to the actual values
#' @description This function plots the predicted values found from running the logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @author Saksham Goel
#' @export
plot_logistic_regression <- function(lr_result)
{
    mean <- vector()

    plot_predictors <- matrix(NA, nrow = 500, ncol = ncol(lr_result$predictors))
    for(i in 1:ncol(lr_result$predictors))
    {
        mean <- append(mean, mean(lr_result$predictors[,i]))
        plot_predictors[, i] <- rep(mean[i], 500)
    }
    for(i in 2:ncol(lr_result$predictors))
    {
        current_predictors <- plot_predictors
        current_predictor <- lr_result$predictors[,i]

        predicted_data <- data.frame(current_predictor = seq(min(current_predictor), max(current_predictor), len = 500))

        current_predictors[, i] <- predicted_data$current_predictor
        predicted_data$response <- 1 / (1 + exp(-(current_predictors %*% lr_result$beta_estimate)))

        #plot(lr_result$response ~ current_predictor)
        #lines(response ~ current_predictor, predicted_data, lwd = 2, col = "green")

        xlabel <- paste("Variable", i)
        plot <- ggplot2::ggplot(x = current_predictor, y = lr_result$response) +
            ggplot2::geom_point(ggplot2::aes(y = lr_result$response, x = current_predictor)) +
            ggplot2::geom_line(ggplot2::aes(y = response, x = current_predictor), data = predicted_data, color = "blue") +
            ggplot2::ylab("Probability") +
            ggplot2::xlab(xlabel)
        print(plot)
    }
}

#' Print Summary
#' @description This function prints summary of the 'logistic_regression' class object created after running the logistic regression.
#' @param lr_result A \code{list} value of class 'logistic_regression', found from running the 'get_beta_estimate' function in this package. The list contains the following objects:
#' \describe{
#'  \item{initial_beta}{The initial values of coefficients.}
#'  \item{beta_estimate}{The estimated coefficients of the logistic regression.}
#'  \item{response}{The \code{double} vector containing the response used for the estimation.}
#'  \item{predictors}{The \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors used for the estimation.}
#' }
#' @author Saksham Goel
#' @export
print_summary <- function(lr_result)
{
    if(class(lr_result) == "logistic_regression")
    {
        cat("Initial values used for beta:", "\n")
        print(lr_result$initial_beta)

        cat("\n", "Estimated coefficients found:", "\n", sep = "")
        print(lr_result$beta_estimate)

        cm <- get_confusion_matrix(lr_result)
        confusion_matrix <- matrix(c(cm$tp, cm$fp, cm$fn, cm$tn), nrow = 2, ncol = 2)
        colnames(confusion_matrix) <- c("Predicted Positive", "Predicted Negative")
        rownames(confusion_matrix) <- c("Actual Positive", "Actual Negative")
        cat("\n", "Confusion Matrix with cut off value = 0.5:", "\n", sep = "")
        print(confusion_matrix)

        cat("\n", "Prevalence: ", get_prevalence(lr_result), "\n", sep = "")
        cat("Accuracy:", get_accuracy(lr_result), "\n")
        cat("Sensitivity:", get_sensitivity(lr_result), "\n")
        cat("Specificity:", get_specificity(lr_result), "\n")
        cat("False Discovery Ratio:", get_false_discovery_ratio(lr_result), "\n")
        cat("Diagnostic Odds Ratio:", get_diagnostic_odds_ratio(lr_result), "\n\n")

        get_beta_confidence_intervals(lr_result)
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
    gender <- sample(c(0, 1), size = n, replace = TRUE)
    age <- round(runif(n, 18, 80))

    predictors <- cbind(gender, age)

    xb <- -9 + 3.5*gender + 0.2*age
    p <- 1/(1 + exp(-xb))
    response <- rbinom(n = n, size = 1, prob = p)

    return(list("response" = response,
                "predictors" = predictors))
}
