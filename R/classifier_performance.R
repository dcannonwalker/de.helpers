# n <- 100
# p <- seq(0,1, length = n)
# tn <- rbinom(n = n, size = 1, prob = p)
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
