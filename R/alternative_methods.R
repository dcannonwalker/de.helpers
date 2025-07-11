#' Fit a standard `edgeR` pipeline for two-group comparison
#' Note to self: to be extended to more general designs
#' @param counts A matrix of counts; columns are samples, rows are tags
#' @param design A design matrix for the columns of `counts`
#' @param ... Additional parameters to pass on
fit_edgeR <- function(counts, design, ...) {
    y <- DGEList(counts = sim$counts,
                 genes = paste0("tag", 1:nrow(sim$counts)),
                 group = trt)
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
}
