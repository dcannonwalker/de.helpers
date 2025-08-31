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
