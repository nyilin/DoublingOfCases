
<!-- badges: start -->

[![DOI](https://img.shields.io/badge/doi-10.1186/s12874.022.01636.3-yellow.svg)](https://doi.org/10.1186/s12874-022-01636-3)
<!-- badges: end -->

DoublingOfCases: Estimating Risk Ratio from Any Standard Epidemiological Design by Doubling the Cases
================

Use the following code to install the package:

``` r
# Package `devtools` required
devtools::install_github("nyilin/DoublingOfCases")
```

This package offers a function `logit_db` that estimates relative risk
in studies of binary outcomes by using the doubling-of-cases approach:

``` r
library(DoublingOfCases)
data("cohort_sim")
# Draw a 1:1 case-control sample from cohort:
i_controls <- sample(which(cohort_sim$y == 0), size = sum(cohort_sim$y),
                     replace = FALSE)
dat_cc <- rbind(
  cbind(cohort_sim[cohort_sim$y == 1, ], w = 1),
  cbind(cohort_sim[i_controls, ], w = sum(cohort_sim$y == 0) / length(i_controls))
)
# Estimate relative risk (RR) by using doubling of cases approach:
# Fit the model:
m_cc <- logit_db(y ~ x + male, data = dat_cc, weight_name = "w")
# Extract model summary:
summ_cc <- summary(m_cc)
# Print the estimates log RR, RR, SEs and p-values:
summ_cc$coefficients
##                    est    exp_est  se_naive se_robust         pval
## (Intercept) -2.3918750 0.09145804 0.1433032 0.1574577 7.730629e-39
## x            0.7016700 2.01711853 0.1841024 0.2078102 8.333284e-04
## male         0.4424209 1.55647068 0.1836683 0.2064312 3.292141e-02
# Compute the 95% CI of the estimated RR:
confint(m_cc, exp = TRUE)
##                   2.5%     97.5%
## (Intercept) 0.06708691 0.1246826
## x           1.34001403 3.0363616
## male        1.03680609 2.3365999
# Extract AIC, log-likelihood and fitted outcomes:
AIC(m_cc)
## [1] 808.3187
logLik(m_cc)
## [1] -401.1593
head(fitted(m_cc))
## [1] 0.1423518 0.2871404 0.2871404 0.2871404 0.2871404 0.2871404
```

**Citation:**

Ning Y, Lam A, Reilly M. Estimating risk ratio from any standard epidemiological design by doubling the cases. BMC Med Res Methodol 22, 157 (2022). https://doi.org/10.1186/s12874-022-01636-3
