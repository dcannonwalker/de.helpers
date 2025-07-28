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
#' @param method One of `c("emp", "emp.paired")`, for two group or
#' two group with paired samples
#' @param n_tags The number of tags to simulate
#' @param ... Arguments passed to sub-functions
#' @export
simulate_effects <- function(n_tags, method, ...) {
    method <- match.arg(method)
    fn <- switch(
        method,
        # dst = simulate_effects.dst,
        # emp = simulate_effects.emp,
        emp.paired = simulated_effects.emp.paired
    )
    effects <- fn(n_tags, ...)
    return(effects)
}

#' Not really supported yet...
#' Simulate effects for a simple two-group design
#' by sampling from estimated coefficients
#' @inheritParams simulate_effects
simulate_effects.emp <- function(n_tags, mean_pars, interval = 1,
                                 min_b1 = log(1.3), ...) {
    if (is.null(mean_pars))
        stop("If method is 'emp', mean_pars must be provided")
    b0 <- sample(mean_pars[, 1], n_tags, replace = TRUE)

    b1 <- sapply(b0, sample_neighborhood, x = b0, y = mean_pars[, 2],
                 interval = interval)
    return(cbind(b0, b1))
}

#' Simulate effects for a two-group design with paired samples
#' by sampling from estimated coefficients
#'
#' Expects that the first column of `mean_pars` is for the intercept,
#' middle columns are for sample effects, and the final column is for
#' treatment effect
#' @inheritParams simulate_effects
#' @inheritParams simulate_counts
#' @param n_pairs The number of sample effects to simulate
#' @param sfx_cols The columns of `mean_pars` to pool for sample effect
#' simulation
#' @param b1_min The smallest absolute value to allow for treatment effects
#' @param p_null The proportion of tags with no treatment effect
simulate_effects.emp.paired <- function(n_tags, n_pairs, mean_pars,
                                        p_null, b1_min = 0,
                                        sfx_cols = NULL, ...) {
    if (is.null(sfx_cols)) {
        message("Pooling all columns except first and last for sample effects")
        sfx_cols <- seq(2, ncol(mean_pars) - 1)
    }
    b1_pars <- mean_pars[, 2]
    n_null <- ceiling(n_tags * p_null)
    b0 <- sample(mean_pars[, 1], n_tags, replace = TRUE)
    b1 <- c(
        rep(0, n_null),
        sample(b1_pars[abs(b1_pars) >= b1_min], n_tags - n_null, replace = TRUE)
    )
    sfx_pars <- c(mean_pars[, sfx_cols])
    sfx <- matrix(sample(sfx_pars, n_tags * n_pairs, replace = TRUE),
                  nrow = n_tags, ncol = n_pairs)
    return(cbind(b0, sfx, b1))
}

#' Not really supported yet...
#' @param de_prob Either a single probability representing the chance
#' that any given tag exhibits differential expression,
#' or a named vector of two probabilities, representing the chance
#' that any given tag exhibits differential expression in the `up` or `down`
#' direction; if a single probability, half of differentially expressed tags
#' will be `up` regulated; used only if `method` is `"dst"`
#' @param theta Exponential parameter to be used if `method` is `"dst"`
#' @param min_log_fc Minimum value added to draws from `"dst"`
simulate_effects.dst <- function(n_tags, theta, de_prob,
                                 distribution = c("exponential"), ...) {
    theta <- theta %||% 1
    de_prob <- de_prob %||% c(up = 0.05, down = 0.05)
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
#' @param mean_pars A matrix of per-tag regression coefficients,
#' with the first column the intercept
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
simulate_counts <- function(mean_pars, dispersion_pars, offset_pars,
                            n_tags, design,
                            offset_options = list(),
                            effects_options = list(),
                            dispersion_interval = log(20)
                            ) {
    n_samples <- nrow(design)
    offset_args <- c(offset_options,
                     list(n_samples = n_samples,
                          n_tags = n_tags,
                          offset_pars = offset_pars))
    offsets <- do.call(simulate_offsets, offset_args)
    effects_args <- c(
        effects_options,
        list(n_tags = n_tags,
             mean_pars = mean_pars)
    )
    effects <- do.call(simulate_effects, effects_args)
    dispersions <- sapply(effects[, 1], sample_neighborhood,
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
            means = means,
            design = design
        )
    )
}

#' Fit `edgeR` to a given dataset, then use the estimated parameters
#' to generate a simulated dataset
#' @param counts A matrix of (presumably real) RNA-Seq counts
#' @param design A design matrix to use to fit `edgeR` model to `counts`
#' @param sim_design A design matrix to use to generate simulated data
#' @inheritParams simulate_counts
#' @param ... Additional arguments
#' to pass to `simulate_counts()`
#' @export
simulate_counts_from_dataset <- function(counts, design,
                                         n_tags,
                                         sim_design,
                                         offset_options,
                                         effects_options, ...) {
    pars <- estimate_parameters(counts = counts, design = design)
    sim <- simulate_counts(mean_pars = pars$mean_pars,
                           dispersion_pars = pars$dispersion_pars,
                           offset_pars = pars$offset_pars, n_tags = n_tags,
                           design = sim_design, offset_options = offset_options,
                           effects_options = effects_options, ...)
    return(
        list(
            pars = pars,
            sim = sim
        )
    )
}
