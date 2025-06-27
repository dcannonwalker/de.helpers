#' Estimate parameters from a real rna-seq data set
#'
#' @param counts a matrix of counts
#' @param group a vector of treatment group labels corresponding to
#' columns of `counts`
#' @param min_count the minimum count to keep a tag in the dataset
#' @param filter_fn function to apply rowwise; result is compared to `min_count`
#' @export
estimate_parameters <- function(counts, group, min_count = 10,
                                filter_fn = sum) {
    idx <- which(apply(counts, 1, filter_fn) > min_count)
    counts <- counts[idx, ]
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    y <- estimateDisp(y)
    fit <- glmFit(y)
    dispersion_pars <- y$tagwise.dispersion
    mean_pars <- fit$coefficients
    list(
        y = y,
        fit = fit,
        offset_pars = y$samples$norm.factors * y$samples$lib.size,
        dispersion_pars = dispersion_pars,
        mean_pars = mean_pars
    )
}
