#' Functions to wrangle `{ngstan}` that I'm not ready to put into
#' the package itself yet
#'
#' Note that `posterior::extract_variable_array()` returns an array that
#' is draws x chains x tags (where applicable) x variable
#'
#' Get the posterior mean for a contrast
#' @param contrast A vector giving the contrast
#' @param comps The `comps` matrix taken from `standata`
#' @param d_pmf The `d_pmf` draws
#' @param beta The `beta` draws
#' @export
get_contrast_posterior_mean <- function(contrast, comps, d_pmf, beta) {
    draws <- fit$draws()
    cp <- ngstan::get_contrast_posterior(contrast, comps, d_pmf, beta)
    return(apply(cp, 1:2, mean))
}
#' Get the posterior mean estimated log2 fold change
#'
#' This currently assumes that there is only a single beta
#' @inheritParams get_contrast_posterior_mean
#' @export
get_log2fc_posterior_mean <- function(comps, d_pfm, beta) {
    # don't drop the dim to keep the same dims as beta
    prob1comp <- d_pmf[, , , which(comps == 1), drop = FALSE]
    log2fc <- log(exp(1), base = 2) * beta * prob1comp
    apply(log2fc, 3, mean)
}
