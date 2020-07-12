#' Double the cases in a dataset
#' @param dat A dataset.
#' @param y_name Name or index of the outcome column in \code{dat}.
#' @return Returns a dataset where each case in the input dataset is doubled. An
#'   ID column (\code{.row_id_db}) is created such that the two rows created
#'   from the same case has the same ID (the name of this new id column begins
#'   with \code{.} so as to minimise the risk of overwriting an existing column
#'   in the input data).
#' @export
double_cases <- function(dat, y_name = "y") {
  dat$.row_id_db <- 1:nrow(dat)
  dat_cases <- dat[dat[, y_name] == 1, ]
  dat_cases[, y_name] <- 0
  dat_db <- rbind(dat_cases, dat)
  dat_db[order(dat_db$.row_id_db), ]
}
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
#' @param weight_name Name or column index in \code{dat} corresponding to the
#'   sampling weights, if \code{dat} is a case-control sample. Default is
#'   \code{NULL}, corresponding to scenarios where the full cohort (or
#'   cross-sectional data) is available.
#' @param dat A dataset (without doubling the cases).
#' @details This function doubles the cases in \code{dat} by using the
#'   \code{\link{double_cases}} function, and then applies a logistic regression
#'   to this modified dataset to estimate the relative risk. The robust SE of
#'   the estimated log relative risk is described in Shouten et al (1993).
#'
#'   If \code{dat} is a (matched) case-control sample, the logistic regression
#'   need to be weighted by the sampling weights corresponding to each subject
#'   in \code{dat}, where the weight is 1 for each case (if all cases are
#'   selected), and the weight for each control is the number of eligible
#'   controls divided by the number of controls sampled (counted within each
#'   sampling stratum if \code{dat} is a matched sample). The robust SE in such
#'   analyses is modified from Shouten et al (1993) by incorporating the
#'   sampling weights in H2.
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
#' @export
logit_db <- function(formula, weight_name = NULL, dat) {
  y_name <- as.character(formula[[2]])
  x_fml <- as.character(formula[[3]])
  if (length(x_fml) > 1) {
    x_fml <- paste("~", paste(x_fml[-1], collapse = x_fml[1]))
  } else {
    x_fml <- paste("~", as.character(x_fml))
  }
  if (!is.null(weight_name)) {
    dat$.weight <- dat[, weight_name]
  } else {
    dat$.weight <- 1
  }
  dat_db <- double_cases(dat = dat, y_name = y_name)
  m <- glm(formula = formula, family = binomial(link = "logit"), data = dat_db,
           weights = .weight)
  coef_mat <- summary(m)$coef
  H1 <- vcov(m)
  # Robust SE
  X <- model.matrix(as.formula(x_fml), data = dat)
  p <- as.vector(exp(X %*% coef(m)))
  y <- unlist(dat[, y_name])
  r <- y - (p / (1 + p)) * (y + 1)
  H2 <- compute_h2(X = X, r = dat$.weight * r)
  vcov_robust <- H1 %*% H2 %*% H1
  se_robust <- sqrt(diag(vcov_robust))
  coef_mat_robust <- data.frame(
    var = rownames(coef_mat),
    est = coef_mat[, 1], se_robust = se_robust,
    pval = 2 * pnorm(abs(coef_mat[, 1] / se_robust), lower.tail = FALSE)
  )
  rownames(coef_mat_robust) <- rownames(coef_mat)
  if (any(p > 1)) {
    warning(simpleWarning(sprintf("%d (%.1f%%) of the %d subjects have estimated probability >1.",
                                  sum(p > 1), mean(p > 1) * 100, nrow(dat))))
    llh <- NA
  } else {
    llh <- sum(dat$.weight * y * log(p) + dat$.weight * (1 - y) * log(1 - p))
  }
  list(obj = m, coef_mat_robust = coef_mat_robust, vcov_robust = vcov_robust,
       fitted_values = as.vector(p), llh = llh,
       AIC = 2 * length(coef(m)) - 2 * llh,
       BIC = log(nrow(dat)) * length(coef(m)) - 2 * llh)
       # gof_pval = generalhoslem::logitgof(obs = y, exp = as.vector(p))$p.value)
}
