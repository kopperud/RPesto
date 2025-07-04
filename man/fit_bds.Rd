% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit_bds}
\alias{fit_bds}
\title{fits the bds model}
\usage{
fit_bds(
  phy,
  sampling_fraction,
  lambda_hat,
  mu_hat,
  eta,
  num_classes = 6,
  sd = 0.587,
  tol = 1e-06,
  condition_survival = TRUE,
  condition_root_speciation = TRUE
)
}
\arguments{
\item{phy}{an object of type phylo}

\item{sampling_fraction}{the probability that each species was sampled in the tree}

\item{lambda_hat}{the overall scale of the log-normal base distribution for the speciation rates. If not specified, the function will estimate it using ML}

\item{mu_hat}{the overall scale of the log-normal base distribution for the extinction rates. If not specified, the function will estimate it using ML}

\item{eta}{the shift rate parameter. If not specified, the function will estimate it using ML}

\item{num_classes}{the number of rate class discretizations (n), such that rate categories is k = n^2}

\item{sd}{the spread parameter for the log-normal base distribution}

\item{tol}{the local error threshold in the numerical ODE solver (per delta_t time step)}

\item{condition_survival}{whether or not to condition on the survival of the left and right lineages descending from the root (default TRUE)}

\item{condition_root_speciation}{whether or not to condition on that there was a speciation event at the root node (default TRUE)}
}
\description{
fits the bds model
}
