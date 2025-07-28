#' Fit a standard `edgeR` pipeline for two-group comparison
#' Note to self: to be extended to more general designs
#' Note that fdr is 'BH' adjustment to p-value
#' @param counts A matrix of counts; columns are samples, rows are tags
#' @param design A design matrix for the columns of `counts`
#' @param ... Additional parameters to pass on
#' @export
fit_edgeR <- function(counts, design, ...) {
    y <- DGEList(counts = counts,
                 genes = paste0("tag", 1:nrow(counts)),
                 group = design[, 2])
    y <- normLibSizes(y)
    y <- estimateDisp(y)
    et <- edgeR::exactTest(y)
    out <- edgeR::topTags(et, n = nrow(et$table), sort.by = "none")
    # these might change when using `glmLRT()` etc.
    expected_colnames <- c(
        "genes",
        "logFC",
        "logCPM",
        "PValue",
        "FDR"
    )
    if (sum(colnames(out) != expected_colnames) != 0)
        warning("topTags colnames not as expected")
    out <- .set_colnames(out)
    return(out)
}

#' Fit a standard `DESeq2` pipeline for two-group comparison
#' Note to self: to be extended to more general designs
#' Note that fdr is 'BH' adjustment to p-value
#' @inheritParams fit_edgeR
#' @export
fit_DESeq2 <- function(counts, design, ...) {
    # second column of design must be treatment group
    trt <- design[, 2]
    y <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = data.frame(
                                            trt = factor(trt)
                                        ), design = ~trt)
    y <- DESeq2::DESeq(y)
    expected_rn <- c("Intercept", "trt_1_vs_0")
    if (DESeq2::resultsNames(y) != expected_rn)
        warning("resultsNames not as expected")
    res <- DESeq2::results(y, name="trt_1_vs_0")
    expected_colnames <- c(
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    )
    if (sum(colnames(res) != expected_colnames) != 0)
        warning("results columns not as expected")
    out <- res
    out[, 'tag'] <- paste0("tag", 1:nrow(counts))
    out <- out[, c(7, 2, 1, 5, 6)]
    out <- .set_colnames(out)
    return(out)
}

#' Fit a standard `limma-trend` or `limma-voom` pipeline for
#' two-group comparison
#' Note to self: to be extended to more general designs
#' Note that fdr is 'BH' adjustment to p-value
#' @inheritParams fit_edgeR
#' @param use_voom Logical - use `limma-voom`? Otherwise use `limma-trend`
#' @export
fit_limma <- function(counts, design, use_voom = TRUE, ...) {
    y <- DGEList(counts = counts,
                 genes = paste0("tag", 1:nrow(sim$counts)),
                 group = design[, 2])
    y <- normLibSizes(y)
    if (use_voom) {
        fit <- .fit_limma(y, design, fn = limma::voom, trend = FALSE)
    } else {
        fit <- .fit_limma(y, design, fn = edgeR::cpm, trend = TRUE,
                          log = TRUE, prior_count = 3)
    }
    out <- limma::topTable(fit, coef = ncol(design))
    expected_colnames <- c(
        "ID",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
    )
    if (sum(colnames(res) != expected_colnames) != 0)
        warning("topTable columns not as expected")
    out <- out[, c(1, 2, 3, 5, 6)]
    out <- .set_colnames(out)
    return(out)
}

.fit_limma <- function(y, design, fn, trend, ...) {
    y_transform <- fn(y, design, ...)
    fit <- limma::lmFit(y_transform, design = design)
    return(limma::eBayes(fit, trend = trend))
}

.set_colnames <- function(out) {
    colnames(out) <- c("tag", "log2fc", "basemean", "prob", "fdr")
    return(out)
}

#' Fit a standard `ngstan` pipeline for two group comparison
#' @inheritParams fit_edgeR
#' @export
fit_ngstan <- function(counts, design, ...) {
    y <- ngstan::seqlist$new(
        counts = counts,
        tags = paste0("tag", 1:nrow(sim$counts))
    )
    y$set_fixed_design(fixed_design = design)
    y$set_mixture_probabilities(c(1, 0.8))
    y$initialize_standata()
    fit <- y$run_model(run_estimation = TRUE, use_multithread = TRUE,
                       grainsize = 125,
                       iter_warmup = 1000,
                       iter_sampling = 1000,
                       parallel_chains = 4,
                       modify_in_place = FALSE)
    draws <- fit$draws()
    comps <- y$standata$comps
    beta <- posterior::extract_variable_array(draws, "beta")
    d_pmf <- posterior::extract_variable_array(draws, "d_pmf")
    logoffset <- apply(
        posterior::extract_variable_array(draws, "log_offset"),
        3, mean)
    prob <- get_contrast_posterior_mean(contrast = 1, comps, d_pmf, beta)
    log2fc <- get_log2fc_posterior_mean(comps, d_pmf, beta)
    out <- data.frame(
        tag = y$tags,
        log2fc = log2fc,
        logoffset = logoffset,
        prob = prob,
        fdr = calc_bfdr(prob)
    )
    return(out)
}
