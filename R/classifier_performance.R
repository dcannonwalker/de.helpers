calc_perf_metrics <- function(p, tn) {
    fpr <- sapply(p, calc_rate, p = p, s = tn)
    tpr <- sapply(p, calc_rate, p = p, s = 1 - tn)
    tfdr <- sapply(p, calc_fdr, p = p, tn = tn)
    cbind(p = p, tn = tn, fpr = fpr, tpr = tpr, tfdr = tfdr)
}
calc_fdr <- function(p0, p, tn) {
    sum(tn[p <= p0]) / sum(p <= p0)
}
calc_rate <- function(p0, p, s) {
    x <- sum(s[p <= p0])
    x / sum(s)
}
#' Calculate the estimated false discovery rate
#' @param p A vector of posterior probabilities
#' @export
calc_bfdr <- function(p) {
    sapply(p, function(p0) {
        sum(p[p <= p0]) / sum(p <= p0)
    })
}

#' Algorithm 1 from Fawcett 2006, for generating ROC points
#' @param p Vector of null probabilities
#' @param tn Vector of actual classes
fawcett1 <- function(p, tn) {
    N <- sum(tn)
    P <- sum(1 - tn)
    ord <- order(p)
    psort <- p[ord]
    tnsort <- tn[ord]
    R <- matrix(nrow = length(p), ncol = 2)
    .fawcett1_internal(R = R, tnsort = tnsort, N = sum(tn), P = P)
}

#' Do the inside step of fawcett1
.fawcett1_internal <- function(R, tnsort, i = 1, FP = 0, TP = 0, N, P) {
    if (i > length(tnsort)) {
        colnames(R) <- c("fpr", "tpr")
        return(R)
    }
    FP <- FP + tnsort[i]
    TP <- TP + (1 - tnsort[i])
    R[i, ] <- c(FP / N, TP / P)
    .fawcett1_internal(R, tnsort, i = i + 1, FP, TP, N, P)
}

#' Algorithm 2 from Fawcett 2006; calculate the AUC
fawcett2 <- function(p, tn) {
    N <- sum(tn)
    P <- sum(1 - tn)
    ord <- order(p)
    tnsort <- tn[ord]
    .fawcett2_internal(tnsort = tnsort, N = N, P = P)
}

#' Do the inside of fawcett2
.fawcett2_internal <- function(A = 0, tnsort, i = 1, FP = 0, TP = 0, N, P) {
    if (i > length(tnsort)) return(A / (N * P))
    A <- A + trapezoid_area(FP, FP + tnsort[i], TP, TP + 1 - tnsort[i])
    .fawcett2_internal(A, tnsort, i = i + 1,
                       FP + tnsort[i], TP + 1 - tnsort[i],
                       N, P)
}

#' Calculate the area of a trapezoidal slice
#' of the AUC
trapezoid_area <- function(x1, x2, y1, y2) {
    abs(x1 - x2) * (y1 + y2) / 2
}

#' Algorithm 3 from Fawcett 2006; vertical averaging of ROC curves
#' @param x0 A vector of fpr values
#' @param rocs A list of ROC curves
fawcett3 <- function(x0, rocs) {
    out <- matrix(nrow = length(x0), ncol = length(rocs) + 1)
    out[, 1] <- x0
    for (i in seq(1, length(x0))) {
        out[i, 2:ncol(out)] <- sapply(rocs, tpr_for_fpr, x0 = x0[i])
    }
    colnames(out) <- c("fpr", paste0("tpr", 1:length(rocs)))
    return(out)
}

#' For a given fpr, find the fpr values on the ROC curve that bracket it
tpr_for_fpr <- function(roc, x0) {
    idx <- suppressWarnings(max(which(roc[, 'fpr'] <= x0)))
    # in case there are 1+ fp before the first tp
    if (idx == -Inf) return(0)
    if (roc[idx, 'fpr'] == x0) return(roc[idx, 'tpr'])
    if (idx == nrow(roc)) return(roc[idx, 'tpr']) # this is a bandaid
    return(
        interpolate(
            x0,
            x1 = roc[idx, 'fpr'],
            x2 = roc[idx + 1, 'fpr'],
            y1 = roc[idx, 'tpr'],
            y2 = roc[idx + 1, 'tpr']
        )
    )
}


#' Interpolate between two points
interpolate <- function(x0, x1, x2, y1, y2) {
    y1 + (x0 - x1) * (y2 - y1) / (x2 - x1)
}
