#' General outline:
#' A generic for simulating counts, paired with S4 classes for
#' each class of input parameters

#' Sample from a neighborhood around a given value
#' @param x0 Given value to compare to `x`
#' @param x Vector from which to get neighborhood indices
#' @param y Vector to sample from
#' @param interval Only pairs where `abs(x - x0) < interval` will
#' be used
#' @export
sample_neighborhood <- function(x0, x, y, interval = 3) {
    distance <- abs(x0 - x)
    yy <- y[distance < interval]
    if (length(yy) == 0) {
        warning("No values within interval")
        return(y[which.min(distance)])
    }
    return(sample(yy, 1))
}

#' Simulate treatment effects
#' @inheritParams simulate_counts
#' @param type One of `c("emp", "dst")`, for sampling from estimated
#' coefficients or sampling from an exponential; if `type` is `"emp"`,
#' `mean_pars` must not be `NULL`
#' @param de_prob Either a single probability representing the chance
#' that any given tag exhibits differential expression,
#' or a named vector of two probabilities, representing the chance
#' that any given tag exhibits differential expression in the `up` or `down`
#' direction; if a single probability, half of differentially expressed tags
#' will be `up` regulated
#' @param theta Exponential parameter to be used if `type` is `"dst"`
#' @param min_log_fc Minimum value added to draws from `"dst"`
#' @param ... Additional arguments passed to sub-functions
simulate_effects <- function(n_tags, mean_pars = NULL, type = c("dst", "emp"),
                             de_prob = c(up = 0.05, down = 0.05),
                             theta = NULL, min_log_fc = log(1.5), ...) {
    type <- match.arg(type)
    if (type == "emp") {
        if (is.null(mean_pars))
            stop("If type is 'emp', mean_pars must be provided")
        effects <- simulate_effects.emp(n_tags = n_tags, mean_pars = mean_pars,
                                        de_prob = de_prob, ...)
    }
    if (type == "dst") {
        theta <- theta %||% 1
        effects <- simulate_effects.dst(n_tags = n_tags, theta = theta,
                                        distribution = "exponential",
                                        de_prob = de_prob,
                                        min_log_fc = min_log_fc, ...)
    }
}

#' Simulate effects by sampling from estimated coefficients
#' @inheritParams simulate_effects
simulate_effects.emp <- function(n_tags, mean_pars, de_prob, interval = 1) {
    # sample base means
    b0 <- sample(mean_pars[, 1], n_tags, replace = TRUE)
    b1 <- sapply(b0, sample_neighborhood, x = b0, y = mean_pars[, 2],
                 interval = interval)

}
#' Simulate RNASeq counts
#' @param mean_pars A matrix of per-tag regression coefficients - should have
#' two columns, with the first column the intercept
#' @param dispersion_pars A vector of estimated per-tag dispersions
#' @param n_tags The number of tags to simulate
#' @param offset_pars A vector of estimated sample offsets
#' @param design The design matrix for the simulated samples - should have
#' two columns, with the first column the intercept
#' @param effects_options A list of options to control the behavior of
#' `simulate_effects()`
simulate_counts <- function(mean_pars, dispersion_pars, n_tags,
                            offset_pars, design) {
    n_samples <- nrow(design)
    effects_args <- c(effects_options, list(n_tags = n_tags))
    if (!is.null(effects_options$type)) {
        if (effects_options$type == "emp")
            effects_args$mean_pars <- mean_pars
    }
    effects <- do.call(simulate_effects, effects_args)
    raw_means <- exp(t(design %*% t(effects)))
}
