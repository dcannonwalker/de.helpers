#' Fit a standard `edgeR` pipeline for two-group comparison
#' Note to self: to be extended to more general designs
#' Note that fdr is 'BH' adjustment to p-value
#' @param counts A matrix of counts; columns are samples, rows are tags
#' @param design A design matrix for the columns of `counts`
#' @param ... Additional parameters to pass on
fit_edgeR <- function(counts, design, ...) {
    y <- DGEList(counts = counts,
                 genes = paste0("tag", 1:nrow(counts)),
                 group = design[, 2])
    y <- edgeR::normLibSizes(y)
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
    if (colnames(out) != expected_colnames)
        warning("topTags colnames not as expected")
    colnames(out) <- c("tag", "log2fc", "log2cpm", "prob", "fdr")
    return(out)
}

#' Fit a standard `DESeq2` pipeline for two-group comparison
#' Note to self: to be extended to more general designs
#' Note that fdr is 'BH' adjustment to p-value
#' @inheritParams fit_edgeR
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
    if (colnames(res) != expected_colnames)
        warning("results columns not as expected")
    out <- res
    out[, 'tag'] <- paste0("tag", 1:nrow(counts))
    out <- out[, c(7, 2, 1, 5, 6)]
    colnames(out) <- c("tag", "log2fc", "basemean", "prob", "fdr")
    return(out)
}
