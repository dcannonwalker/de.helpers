#' Construct an ROC curve for a given model fit
#' @param method_data Model fit for a single method,
#' as read in by `read_method_data()`
#' @export
make_roc_curve <- function(method_data) {
    fawcett1(p = method_data[["prob"]],
             tn = method_data[["true_null"]])
}

#' Construct the average ROC curve for a set of model fits
#' @param method_data List of model fits for a single method
#' @param x0 A vector of fpr or nominal fdr values at which to take
#' the averages of the ROC or FDR curves
#' @export
make_average_roc_curve <- function(method_data, x0) {
    roc_curves <- lapply(method_data, make_roc_curve)
    fawcett3(x0, roc_curves)
}

#' Construct an FDR curve for a given model fit
#' @inheritParams make_roc_curve
#' @param misleading_names Output will have misleading column names
#' so that it can be passed to `fawcett3()` without having to change
#' the column names later (or that function to be more generic);
#' change it later if you feel like it
#' @export
make_fdr_curve <- function(method_data, misleading_names = FALSE) {
    p <- method_data[["prob"]]
    fdr <- method_data[["fdr"]]
    tfdr <- sapply(p, calc_fdr, p = p, tn = method_data[["true_null"]])
    ord <- order(p)
    if (misleading_names) {
        cbind(fpr = fdr[ord], tpr = tfdr[ord])
    } else cbind(fdr = fdr[ord], tfdr = tfdr[ord])
}

#' Construct the average FDR curve for a set of model fits
#' @inheritParams make_average_roc_curve
#' @export
make_average_fdr_curve <- function(method_data, x0) {
    fdr_curves <- lapply(method_data, make_fdr_curve, misleading_names = TRUE)
    average_fdr <- fawcett3(x0, fdr_curves)
    colnames(average_fdr) <- c("fdr", paste0("tfdr", 1:length(rocs)))
    average_fdr
}
