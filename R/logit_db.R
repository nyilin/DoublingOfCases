#' Double the cases in a dataset
#' @param data A dataset.
#' @param y_name Name or index of the outcome column in \code{data}.
#' @return Returns a dataset where each case in the input dataset is doubled. An
#'   ID column (\code{.row_id_db}) is created such that the two rows created
#'   from the same case has the same ID (the name of this new id column begins
#'   with \code{.} so as to minimise the risk of overwriting an existing column
#'   in the input data).
double_cases <- function(data, y_name = "y") {
  data$.row_id_db <- 1:nrow(data)
  dat_cases <- data[data[, y_name] == 1, ]
  dat_cases[, y_name] <- 0
  dat_db <- rbind(dat_cases, data)
  dat_db[order(dat_db$.row_id_db), ]
}
#' Compute H2 matrix for robust covariance
compute_h2 <- function(X, r) {
  nr <- nrow(X)
  nc <- ncol(X)
  h2 <- do.call("rbind", lapply(1:nr, function(i) {
    h2_mat <- matrix(X[i, ], ncol = 1) %*% matrix(X[i, ], nrow = 1) * (r[i] ^ 2)
    as.vector(h2_mat)
  }))
  matrix(apply(h2, 2, sum), ncol = nc)
}
#' Estimate the relative risk from a binary outcome by using the
#' doubling-of-cases approach
#' @param formula Formula.
#' @param weight_name Name or column index in \code{data} corresponding to the
#'   sampling weights, if \code{data} is a case-control sample. Default is
#'   \code{NULL}, corresponding to scenarios where the full cohort (or
#'   cross-sectional data) is available.
#' @param data A \code{data.frame} or named matrix.
#' @param ... Other input to \code{\link{glm}}.
#' @details This function doubles the cases in \code{data}, and then applies a
#'   logistic regression to this modified dataset to estimate the relative risk.
#'   The robust SE of the estimated log relative risk is described in Schouten et
#'   al (1993).
#'
#'   If \code{data} is a (matched) case-control sample, the logistic regression
#'   need to be weighted by the sampling weights corresponding to each subject
#'   in \code{data}, where the weight is 1 for each case (if all cases are
#'   selected), and the weight for each control is the number of eligible
#'   controls divided by the number of controls sampled (counted within each
#'   sampling stratum if \code{data} is a matched sample). The robust SE in such
#'   analyses is modified from Schouten et al (1993) by incorporating the
#'   sampling weights in H2.
#' @return Returns an object of class \code{db} that inherits from the
#'   \code{\link{glm}} class. Note that although \code{family} of this object is
#'   logistic, the underlying model fitted is log-binomial. Items \code{call},
#'   \code{deviance}, \code{null.deviance} and \code{weights} are not available.
#' @references
#' \itemize{
#'   \item{Schouten, E.G., Dekker, J.M., Kok, F.J., Cessie, S.L., Van
#'   Houwelingen, H.C., Pool, J. and Vandnbroucke, J.P. (1993). Risk ratio and
#'   rate ratio estimation in case‐cohort designs: Hypertension and
#'   cardiovascular mortality. Statist. Med., 12: 1733-1745.
#'   doi:10.1002/sim.4780121808}
#'   \item{Knol, M. J., Le Cessie, S., Algra, A., Vandenbroucke, J. P., &
#'   Groenwold, R. H. (2012). Overestimation of risk ratios by odds ratios in
#'   trials and cohort studies: alternatives to logistic regression. CMAJ :
#'   Canadian Medical Association journal, 184(8), 895–899.
#'   https://doi.org/10.1503/cmaj.101715}
#' }
#' @examples
#' library(DoublingOfCases)
#' data("cohort_sim")
#' # Draw a 1:1 case-control sample from cohort:
#' i_controls <- sample(which(cohort_sim$y == 0), size = sum(cohort_sim$y),
#'                      replace = FALSE)
#' dat_cc <- rbind(
#'   cbind(cohort_sim[cohort_sim$y == 1, ], w = 1),
#'   cbind(cohort_sim[i_controls, ], w = sum(cohort_sim$y == 0) / length(i_controls))
#' )
#' # Estimate relative risk (RR) by using doubling of cases approach:
#' # Fit the model:
#' m_cc <- logit_db(y ~ x + male, data = dat_cc, weight_name = "w")
#' # Extract model summary:
#' summ_cc <- summary(m_cc)
#' # Print the estimates log RR, RR, SEs and p-values:
#' summ_cc$coefficients
#' # Compute the 95% CI of the estimated RR:
#' confint(m_cc, exp = TRUE)
#' # Extract AIC, log-likelihood and fitted outcomes:
#' AIC(m_cc)
#' logLik(m_cc)
#' fitted(m_cc)
#' @export
logit_db <- function(formula, weight_name = NULL, data, ...) {
  y_name <- as.character(formula[[2]])
  x_fml <- as.character(formula[[3]])
  if (length(x_fml) > 1) {
    x_fml <- paste("~", paste(x_fml[-1], collapse = x_fml[1]))
  } else {
    x_fml <- paste("~", as.character(x_fml))
  }
  if (!is.null(weight_name)) {
    data$.weight <- data[, weight_name]
  } else {
    data$.weight <- 1
  }
  dat_db <- double_cases(data = data, y_name = y_name)
  m <- glm(formula = formula, family = binomial(link = "logit"), data = dat_db,
           weights = .weight, ...)
  # Get fitted values and residuals of log-binomial
  y <- unlist(data[, y_name])
  X <- model.matrix(as.formula(x_fml), data = data)
  m$y <- y
  m$x <- X
  m$model <- m$model[!duplicated(dat_db$.row_id_db), ]
  p <- as.vector(exp(X %*% coef(m)))
  m$residuals <- y - p
  m$linear.predictors <- log(p)
  m$fitted.values <- p
  m$df.null <- nrow(data)
  m$df.residual <- nrow(data) - length(coef(m))
  m$prior.weights <- data$.weight
  if (any(p > 1)) {
    warning(simpleWarning(sprintf(
      "%d (%.1f%%) of the %d subjects have estimated probability >1.",
      sum(p > 1), mean(p > 1) * 100, nrow(data)
    )))
    llh <- NA
  } else {
    llh <- sum(data$.weight * y * log(p) + data$.weight * (1 - y) * log(1 - p))
  }
  m$log_likelihood <- llh
  m$aic <- 2 * length(coef(m)) - 2 * llh
  m$H1 <- vcov(m)
  # The following are not yet implemented:
  m$call <- NULL #**** Not sure if can make this the correct call to logit_db
  m$deviance <- NULL
  m$null.deviance
  m$weights <- NULL
  class(m) <- c("db", "glm", "lm")
  m
}
#' @describeIn summary.db
#' Extract corrected covariance matrix
#' @export
vcov.db <- function(object, ...) {
  r <- object$y - (object$fitted.values / (1 + object$fitted.values)) * (object$y + 1)
  H2 <- compute_h2(X = object$x, r = object$prior.weights * r)
  object$H1 %*% H2 %*% object$H1
}
#' @describeIn summary.db
#' Extract log-likelihood with respect to the log-binomial model of interest
#' @export
logLik.db <- function(object, ...) {
  object$log_likelihood
}
#' @describeIn summary.db
#' Extract AIC
#' @export
AIC.db <- function(object, ..., k = 2) {
  k * length(coef(object)) - 2 * logLik(object)
}
#' @describeIn summary.db
#' Confidence intervals for estimated coefficients
#' @export
confint.db <- function(object, parm, exp = FALSE, level = 0.95, ...) {
  level <- as.numeric(level)
  if (is.na(level) || (level <= 0 | level >= 1)) {
    stop(simpleError("level should be between 0 and 1."))
  }
  vcov_robust <- vcov.db(object)
  se_robust <- sqrt(diag(vcov_robust))
  est <- coef(object)
  if (!missing(parm)) {
    est <- est[parm]
    se_robust <- se_robust[parm]
  }
  p <- (1 - level) / 2
  q <- qt(p = p, df = object$df.residual)
  ci_mat <- data.frame(l = est - q * se_robust,
                       u = est + q * se_robust)
  names(ci_mat) <- c(paste0(round(p * 100, 2), "%"),
                     paste0(round((1 - p) * 100, 2), "%"))
  rownames(ci_mat) <- names(est)
  if (exp) {
    exp(ci_mat)
  } else {
    ci_mat
  }
}
#' Summarise the fitted model for the doubling of cases approach
#' @param object The model fitted using \code{\link{logit_db}}.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param exp Whether to generate confidence interval for the exponentiated
#'   coefficients (i.e., the estimated risk ratio, if \code{TRUE}, the default),
#'   or for the estimated coefficients (if \code{FALSE}).
#' @param level The confidence level required. Default is 0.95 (i.e., 95\%).
#' @param ... Additional input.
#' @details See usage in \code{\link{logit_db}}.
#' @export
summary.db <- function(object, ...) {
  vcov_robust <- vcov.db(object)
  se_robust <- sqrt(diag(vcov_robust))
  coefficients <- data.frame(
    est = coef(object), exp_est = exp(coef(object)),
    se_naive = sqrt(diag(object$H1)), se_robust = se_robust,
    pval = 2 * pt(q = abs(coef(object) / se_robust), df = object$df.residual,
                  lower.tail = FALSE)
  )
  list(call = object$call, aic = object$aic,
       df = c(object$rank, object$df.residual, length(coef(object))),
       coefficients = coefficients,
       cov.unscaled = vcov_robust, cov.scaled = vcov_robust)
}
