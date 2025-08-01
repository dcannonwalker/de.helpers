% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_counts.R
\name{simulate_effects.emp.paired}
\alias{simulate_effects.emp.paired}
\title{Simulate effects for a two-group design with paired samples
by sampling from estimated coefficients}
\usage{
simulate_effects.emp.paired(
  n_tags,
  n_pairs,
  mean_pars,
  p_null,
  b1_min = 0,
  sfx_cols = NULL,
  ...
)
}
\arguments{
\item{n_tags}{The number of tags to simulate}

\item{n_pairs}{The number of sample effects to simulate}

\item{mean_pars}{A matrix of per-tag regression coefficients,
with the first column the intercept}

\item{p_null}{The proportion of tags with no treatment effect}

\item{b1_min}{The smallest absolute value to allow for treatment effects}

\item{sfx_cols}{The columns of \code{mean_pars} to pool for sample effect
simulation}

\item{...}{Arguments passed to sub-functions}
}
\description{
Expects that the first column of \code{mean_pars} is for the intercept,
middle columns are for sample effects, and the final column is for
treatment effect
}
