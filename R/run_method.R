#' Fit a model to saved simulated data
#' @inheritParams make_sim_dir
#' @inheritParams make_dataset_dir
#' @param method The differential expression analysis method to run
#' @param ... Additional arguments passed to the selected method
#' @export
run_method <- function(simulation_id, dataset_id, method = c("edgeR", "DESeq2", "limma", "ngstan"),
                       root = file.path("out", "simulation_studies"), ...) {
    method <- match.arg(method)
    design <-
        as.matrix(
            read.table(file.path(root, simulation_id, "design"), header = TRUE)
        )
    counts <- read.table(file.path(root, simulation_id, dataset_id, "counts"), header = TRUE)
    method_fn <- switch(method,
                        edgeR = fit_edgeR,
                        DESeq2 = fit_DESeq2,
                        limma = fit_limma,
                        ngstan = fit_ngstan)
    fit <- method_fn(counts = counts, design = design, ...)
    message(
        glue::glue("Saving fit to",
                   " {file.path(root, simulation_id, dataset_id, paste0(method, '.qs2'))}")
    )
    if (names(fit) == c("fit", "fit_raw")) {
        qs2::qs_save(fit$fit, file.path(root, simulation_id, dataset_id, paste0(method, ".qs2")))
        qs2::qs_save(fit$fit_raw, file.path(root, simulation_id, dataset_id, paste0(method, "_raw.qs2")))
    } else {
        qs2::qs_save(fit, file.path(root, simulation_id, dataset_id, paste0(method, ".qs2")))
    }
}
