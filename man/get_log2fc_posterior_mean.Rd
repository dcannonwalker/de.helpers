% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngstan.R
\name{get_log2fc_posterior_mean}
\alias{get_log2fc_posterior_mean}
\title{Get the posterior mean estimated log2 fold change}
\usage{
get_log2fc_posterior_mean(comps, d_pfm, beta)
}
\arguments{
\item{comps}{The \code{comps} matrix taken from \code{standata}}

\item{beta}{The \code{beta} draws}
}
\description{
This currently assumes that there is only a single beta
}
