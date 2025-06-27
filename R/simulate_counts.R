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
#' @param method One of `c("emp", "dst")`, for sampling from estimated
#' coefficients or sampling from an exponential; if `method` is `"emp"`,
#' `mean_pars` must not be `NULL`
#' @param de_prob Either a single probability representing the chance
#' that any given tag exhibits differential expression,
#' or a named vector of two probabilities, representing the chance
#' that any given tag exhibits differential expression in the `up` or `down`
#' direction; if a single probability, half of differentially expressed tags
#' will be `up` regulated; used only if `method` is `"dst"`
#' @param theta Exponential parameter to be used if `method` is `"dst"`
#' @param min_log_fc Minimum value added to draws from `"dst"`
#' @param ... Additional arguments passed to sub-functions
#' @export
simulate_effects <- function(n_tags, mean_pars = NULL, method = c("dst", "emp"),
                             de_prob = NULL,
                             theta = NULL, min_log_fc = log(1.5), ...) {
    method <- match.arg(method)
    if (method == "emp") {
        if (is.null(mean_pars))
            stop("If method is 'emp', mean_pars must be provided")
        effects <- simulate_effects.emp(n_tags = n_tags, mean_pars = mean_pars,
                                        ...)
    }
    if (method == "dst") {
        theta <- theta %||% 1
        de_prob <- de_prob %||% c(up = 0.05, down = 0.05)
        b1 <- simulate_effects.dst(n_tags = n_tags, theta = theta,
                                   distribution = "exponential",
                                   de_prob = de_prob,
                                   min_log_fc = min_log_fc, ...)
    }
    return(effects)
}

#' Simulate effects by sampling from estimated coefficients
#' @inheritParams simulate_effects
simulate_effects.emp <- function(n_tags, mean_pars, interval = 1) {
    b0 <- sample(mean_pars[, 1], n_tags, replace = TRUE)

    b1 <- sapply(b0, sample_neighborhood, x = b0, y = mean_pars[, 2],
                 interval = interval)
    return(cbind(b0, b1))
}

simulate_effects.dst <- function(n_tags, theta, de_prob,
                                 distribution = c("exponential")) {
    distribution <- match.arg(distribution)
    if (sum(de_prob) < 0 || sum(de_prob) > 1)
        stop("sum(de_prob) must be between 0 and 1")
    if (length(de_prob) == 1) {
        message("simulating equal chance of up and down regulation")
        de_prob <- c(de_prob / 2, de_prob / 2)
    }
    de <- sample(c(-1, 0, 1), size = n_tags, replace = TRUE,
                 prob = c(de_prob['down'], 1 - sum(de_prob), de_prob['up']))
}

#' Simulate per-sample offsets
#' @inheritParams simulate_effects
#' @inheritParams simulate_counts
#' @param n_samples The number of samples for which to simulate offsets
simulate_offsets <- function(offset_pars, n_samples, n_tags, method = c("default")) {
    method <- match.arg(method)
    edgeR::makeCompressedMatrix(
        sample(offset_pars, n_samples, replace = TRUE),
        dims = c(n_tags, n_samples)
    )
}

#' Simulate RNASeq counts
#' @param mean_pars A matrix of per-tag regression coefficients - should have
#' two columns, with the first column the intercept
#' @param dispersion_pars A vector of estimated per-tag dispersions
#' @param dispersion_interval After simulating the base mean `b0`,
#' the dispersion is sampled from
#' `dispersion_pars[abs(mean_pars[, 1] - b0) < dispersion_interval]`
#' @param n_tags The number of tags to simulate
#' @param design The design matrix for the simulated samples - should have
#' two columns, with the first column the intercept
#' @param offset_pars A vector of estimated sample offsets
#' @param offset_options A list of options to control the behavior of
#' `simulate_offsets()`
#' @param effects_options A list of options to control the behavior of
#' `simulate_effects()`
#' @export
simulate_counts <- function(mean_pars, dispersion_pars,
                            dispersion_interval = log(20), n_tags,
                            design,
                            offset_pars,
                            offset_options = list(),
                            effects_options = list()) {
    n_samples <- nrow(design)
    offset_args <- c(offset_options,
                     list(n_samples = n_samples,
                          n_tags = n_tags,
                          offset_pars = offset_pars))
    offsets <- do.call(simulate_offsets, offset_args)
    effects_args <- c(effects_options, list(n_tags = n_tags))
    if (!is.null(effects_options$method)) {
        if (effects_options$method == "emp")
            effects_args$mean_pars <- mean_pars
    }
    effects <- do.call(simulate_effects, effects_args)
    dispersions <- sapply(mean_pars[, 1], sample_neighborhood,
                          x = mean_pars[, 1], y = dispersion_pars,
                          interval = dispersion_interval)
    means <- exp(offsets + t(design %*% t(effects)))
    counts <- matrix(nrow = n_tags, ncol = n_samples)
    for (s in seq(1, n_samples)) {
        counts[, s] <- rnbinom(n_tags, mu = means[, s], size = 1 / dispersions)
    }
    return(
        list(
            counts = counts,
            effects = effects,
            offsets = offsets,
            dispersions = dispersions,
            means = means
        )
    )
}
