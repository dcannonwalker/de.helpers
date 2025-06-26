#' General outline:
#' A generic for simulating counts, paired with S4 classes for
#' each class of input parameters

#' Sample a dispersion given a mean and a set of mean and dispersion pairs
#' @param x Given mean
#' @param mean_pars Vector of estimated means
#' @param dispersion_pars Vector of estimated dispersions
#' @param interval Only pairs where `abs(mean_pars - x) < interval` will
#' be used
#' @export
sample_dispersion <- function(x, mean_pars, dispersion_pars, interval = 20) {
    distance <- abs(x - mean_pars)
    y <- dispersion_pars[distance < interval]
    if (length(y) == 0) {
        warning("No means within interval")
        return(dispersion_pars[which.min(distance)])
    }
}
