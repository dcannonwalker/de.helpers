% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_counts.R
\name{simulate_counts_from_dataset}
\alias{simulate_counts_from_dataset}
\title{Fit \code{edgeR} to a given dataset, then use the estimated parameters
to generate a simulated dataset}
\usage{
simulate_counts_from_dataset(
  counts,
  design,
  n_tags,
  sim_design,
  offset_options,
  effects_options,
  full_output = FALSE,
  ...
)
}
\arguments{
\item{counts}{A matrix of (presumably real) RNA-Seq counts}

\item{design}{A design matrix to use to fit \code{edgeR} model to \code{counts}}

\item{n_tags}{The number of tags to simulate}

\item{sim_design}{A design matrix to use to generate simulated data}

\item{offset_options}{A list of options to control the behavior of
\code{simulate_offsets()}}

\item{effects_options}{A list of options to control the behavior of
\code{simulate_effects()}}

\item{...}{Additional arguments
to pass to \code{simulate_counts()}}
}
\description{
Fit \code{edgeR} to a given dataset, then use the estimated parameters
to generate a simulated dataset
}
