#' Estimate parameters from a real rna-seq data set
#'
#' @param counts
#' @param group
estimate_parameters <- function(counts, group, min_count = 10,
                                filter_fn = sum, ...) {
    idx <- which(apply(counts, 1, filter_fn) < min_count)
    counts <- counts[idx, ]
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    y <- estimateDisp(y)
    dispersion_pars <- y$tagwise.dispersion
    mean_pars <- apply(counts, 1, mean)
    list(
        dispersion_pars = dispersion_pars,
        mean_pars = mean_pars
    )
}
