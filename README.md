
<!-- badges: start -->

[![DOI](https://img.shields.io/badge/doi-10.1186/s12874.022.01636.3-yellow.svg)](https://doi.org/10.1186/s12874-022-01636-3)
<!-- badges: end -->

DoublingOfCases: Estimating Risk Ratio from Any Standard Epidemiological
Design by Doubling the Cases
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
## (Intercept) -2.3480440 0.09555588 0.1443916 0.1595498 4.378952e-37
## x            0.8401388 2.31668852 0.1834710 0.2066872 6.181787e-05
## male         0.2404291 1.27179470 0.1830716 0.2053742 2.426756e-01
# Compute the 95% CI of the estimated RR:
confint(m_cc, exp = TRUE)
##                   2.5%     97.5%
## (Intercept) 0.06980476 0.1308066
## x           1.54242984 3.4796044
## male        0.84894015 1.9052718
# Extract AIC, log-likelihood and fitted outcomes:
AIC(m_cc)
## [1] 807.852
logLik(m_cc)
## [1] -400.926
head(fitted(m_cc))
## [1] 0.1215275 0.2815413 0.2815413 0.2815413 0.2815413 0.2815413
```

**Citation:**

Ning Y, Lam A, Reilly M. Estimating risk ratio from any standard
epidemiological design by doubling the cases. BMC Med Res Methodol 22,
157 (2022). <https://doi.org/10.1186/s12874-022-01636-3>
