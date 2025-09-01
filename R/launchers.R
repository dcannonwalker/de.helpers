#' Prepare a system prompt
#'
#' @export
prepare_slurm_prompt <- function(simulation_id, dataset_id, method, slurm_script_name, ..., root = "out/simulation_studies") {
    args <- c(...)
    job_name <- glue::glue("{simulation_id}_{dataset_id}_{method}")
    system_args <- c(glue::glue("-o {file.path(root, simulation_id, 'logs', job_name)}"),
                     glue::glue("--job-name {job_name}"),
                     "slurm/single_method_launcher",
                     simulation_id,
                     dataset_id,
                     method)
    if (length(args) >= 1) {
        system_args <- c(system_args, args)
    }
    system_args
}

#' Launch single method using system prompt
#'
#' For `ngstan`, the first two elements of `...` will be used
#' as `iter_warmup` and `iter_sampling`
#' @export
launch_single_method <- function(simulation_id, dataset_id, method, ..., root = "out/simulation_studies") {
    system_args <- prepare_slurm_prompt(simulation_id, dataset_id, method,
                                        slurm_script_name = "slurm/single_method_launcher",
                                        ..., root = root)
    system2("sbatch", args = system_args)
}

#' Launch model fits for a set of methods for all data sets
#' in a given simulation study
#' @export
launch_simulation_study <- function(simulation_id,
                                    methods = c(
                                        "ngstan",
                                        "edgeR",
                                        "DESeq2",
                                        "limma"
                                    ), ...,
                                    verbose = TRUE,
                                    root = file.path("out", "simulation_studies")) {
    args <- c(...)
    methods <- match.arg(methods, several.ok = TRUE)
    dataset_ids <- read.table(file.path(root, simulation_id, "dataset_ids.txt"), header = TRUE)
    launcher_args <- list(
        simulation_id = simulation_id,
        root = root
    )
    for (id in dataset_ids$dataset_id) {
        launcher_args[["dataset_id"]] <- id
        for (method in methods) {
            launcher_args[["method"]] <- method
            if (method == "ngstan" && length(args) == 2) {
                launcher_args[["iter_warmup"]] <- args[1]
                launcher_args[["iter_sampling"]] <- args[2]
            }
            if (verbose) message(glue::glue("Launching {method} for data set {id}..."))
            do.call(launch_single_method, args = launcher_args)
        }
    }
}
