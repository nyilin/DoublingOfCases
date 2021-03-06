% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logit_db.R
\name{compute_h2}
\alias{compute_h2}
\title{Compute H2 matrix for robust covariance}
\usage{
compute_h2(X, r)
}
\arguments{
\item{X}{Design matrix from original data (without doubling the cases).}

\item{r}{(Weighted) residuals for each subject.}
}
\description{
Compute H2 matrix for robust covariance
}
