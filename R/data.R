#' Simulated cohort
#'
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#'   \item{y}{A binary outcome.}
#'   \item{x}{A binary exposure.}
#'   \item{male}{Gender of subjects (1 for male and 0 for female).}
#' }
#'
#' In this simulated cohort, the probability of being male (z=1) was 0.4. The
#' probability of being exposed (x=1) was 0.4 for male and 0.2 for female. The
#' outcome was generated from the following log-binomial model:
#' log{P(Y=1)} = log(0.1) + log(2)*x + log(1.5)*z.
#'
"cohort_sim"
