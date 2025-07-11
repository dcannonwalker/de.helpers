#' Functions to wrangle `{ngstan}` that I'm not ready to put into
#' the package itself yet
#'
#' Get the posterior mean for a contrast
#' @param contrast A vector giving the contrast
#' @param comps The `comps` matrix taken from `standata`
#' @param d_pmf The `d_pmf` draws
#' @param beta The `beta` draws
get_contrast_posterior_mean <- function(contrast, comps, d_pmf, beta) {
    draws <- fit$draws()
    cp <- ngstan::get_contrast_posterior(contrast, comps, d_pmf, beta)
    return(apply(cp, 1:2, mean))
}
#' Get the posterior mean estimated log2 fold change
#'
#' This currently assumes that there is only a single
get_log2fc_posterior_mean <- function(d_pfm, beta) {

}
