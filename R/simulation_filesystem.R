#' Solicit yes/no answer from user
#' @param x Message with the question
yesno <- function(x) {
    cli::cli_inform(x)
    out <- utils::menu(c("yes", "no"))
    out  %in% c("yes", 1)
}

#' Set up a directory to contain all files for a given simulation study
#' @param simulation_id The simulation study identifier
#' @param root The parent directory for all simulation studies
#' @export
make_sim_dir <- function(simulation_id,
                         root = file.path("out", "simulation_studies")) {
    path <- file.path(root, simulation_id)
    exists <- dir.exists(path)
    if (exists) {
        message("Directory already exists; not creating")
    } else {
        msg <- glue::glue("Creating directory at",
                          " {path};",
                          " do you wish to proceed?")
        sel <- yesno(msg)
        if (sel) dir.create(path, recursive = TRUE)
    }
    doc <- yesno("Add a line of metadata to simulation study index?")
    if (doc) {
        docfile <- file.path(root, "study_index.txt")
        docexists <- file.exists(docfile)
        append <- FALSE
        col.names <- TRUE
        if (docexists) {
            append <- TRUE
            col.names <- FALSE
        }
        desc <- readline("Enter a one-line description: ")
        row <- data.frame(simulation_id = simulation_id, description = desc)
        message(glue::glue("Writing description to {docfile}"))
        write.table(row, file = docfile, append = append,
                    col.names = col.names, row.names = FALSE)
    }
    return(path)
}

#' Set up a directory to contain all the files for a single simulated data set
#' within a simulation study
#' @param dataset_id The data set identifier
#' @inheritParams make_sim_dir
#' @export
make_dataset_dir <- function(dataset_id, simulation_id,
                             root = file.path("out", "simulation_studies")) {
    path <- file.path(root, simulation_id, dataset_id)
    exists <- dir.exists(path)
    if (exists) {
        message("Directory already exists; not creating")
    } else {
        if (interactive()) {
            msg <- glue::glue("Creating directory at",
                              " {path};",
                              " do you wish to proceed?")
            sel <- yesno(msg)
        } else sel <- TRUE
        if (sel) dir.create(path, recursive = TRUE)
    }
    docfile <- file.path(root, simulation_id, "dataset_ids.txt")
    docexists <- file.exists(docfile)
    append <- FALSE
    col.names <- TRUE
    if (docexists) {
        append <- TRUE
        col.names <- FALSE
    }
    row <- data.frame(dataset_id = dataset_id)
    write.table(row, file = docfile, append = append,
                col.names = col.names, row.names = FALSE)
    return(path)
}

#' Generate an identifier
#' @param n_ids Number of ids to generate
#' @param n_char Number of characters
#' @param ids_to_check A vector of other ids to check against
#' @param seed Seed value to set
#' @export
generate_id <- function(n_ids = 1, n_char = 10, ids_to_check = NULL,
                        seed = NULL, i = 0) {
    if (i > 5) stop("i > 5; are you stuck in a loop?")
    if (n_ids > 1) {
        ids <- sapply(1:n_ids, function(i) {
            generate_id(n_char = n_char, ids_to_check = ids_to_check,
                        seed = NULL, i = 0)
        })
        if (length(unique(ids)) < length(ids)) {
            warning("non-unique ids generated")
        }
        return(ids)
    }
    if (!is.null(seed)) set.seed(seed)
    id <- paste0(sample(letters, n_char, replace = TRUE), collapse = "")
    if (!is.null(ids_to_check)) {
        exists <- id %in% ids_to_check
        if (exists) {
            message("id already in use; regenerating")
            id <- generate_id(n_char = n_char, ids_to_check = ids_to_check,
                              i = i + 1)
        }
    }
    return(id)
}

#' Save simulated data sets for a simulation study
#' @inheritParams make_sim_dir
#' @param sim_data A list of simulated data sets, each element named
#' matching one of the data set ids in the simulation study
#' @param design The design matrix for the simulated counts
#' @export
save_datasets <- function(simulation_id, sim_data, design,
                          root = file.path("out", "simulation_studies")) {
    simulation_root <- file.path(root, simulation_id)
    message(
        glue::glue("Saving design matrix to",
                   " {file.path(simulation_root, 'design')}...")
    )
    write.table(design, file = file.path(simulation_root, "design"))
    lapply(names(sim_data), function(dataset_id) {
        dataset_root <- file.path(simulation_root, dataset_id)
        lapply(names(sim_data[[dataset_id]]), function(tbl) {
            message(
                glue::glue("Saving {tbl} to",
                           " {file.path(dataset_root, tbl)}...")
            )
            write.table(sim_data[[dataset_id]][[tbl]],
                        file.path(dataset_root, tbl))
        })
    })
    message("Simulated data sets saved")
}

#' Read in all the saved fits for a given simulation study and method
#' @inheritParams make_sim_dir
#' @inheritParams make_dataset_dir
#' @param method The name of the method whose fits we want to retrieve
#' @export
read_method_data <- function(
        simulation_id,
        root = "out/simulation_studies",
        method = c("edgeR",
                   "DESeq2",
                   "limma",
                   "ngstan")
) {
    method <- match.arg(method)
    dataset_ids <- read.table(
        file.path(root, simulation_id, "dataset_ids.txt"),
        header = TRUE
    )
    out_list <- lapply(dataset_ids$dataset_id, .read_method_data,
                       simulation_id = simulation_id, method = method, root = root)

}

.read_method_data <- function(dataset_id, simulation_id, root, method) {
    method_data <- qs2::qs_read(
        file.path(root, simulation_id, dataset_id,
                  glue::glue("{method}.qs2"))
    )
    effects <- read.table(
        file.path(root, simulation_id, dataset_id,
                  "effects"), row.names = NULL
    )
    true_null <- effects[, "b1"] == 0
    data.frame(simulation_id = simulation_id, dataset_id = dataset_id,
               true_null = true_null, method_data)
}

#' Get the standard file path for saved ROC or FDR curve data
#' @inheritParams save_method_curves
get_curve_path <- function(simulation_id, root, method, type) {
    file.path(root, simulation_id, glue::glue("{method}_{type}_curve"))
}

#' A wrapper to create and save average ROC and FDR curves
#' @section Running and analyzing simulation studies:
#' To use the helpers to run and analyze a simulation study, you're intended
#' to follow these steps in order:
#'
#' 1. Simulate data sets and save them
#' 2. Fit the set of methods and save their outputs for each data set
#' 3. Create tables for average ROC and FDR curves and save them
#' 4. Plot the ROC and FDR curves and save the plots
#' @section Simulate data:
#' `generate_id()`, `make_sim_dir()`, `make_dataset_dir()`, `simulate_counts_from_dataset()`,
#' `save_datasets()`
#' @section Fit methods:
#' `run_method()`
#' @section Create summary tables:
#' `save_method_curves()`
#' @section Save plots:
#' `save_curve_plots()`
#' @inheritParams read_method_data
#' @inheritParams make_average_roc_curve
#' @param type ROC, FDR, or both?
#' @export
save_method_curves <- function(
        simulation_id,
        root = "out/simulation_studies",
        method = c("edgeR",
                   "DESeq2",
                   "limma",
                   "ngstan"),
        type = c("roc", "fdr"),
        x0 = seq(0, 1, length = 100)
) {
    method <- match.arg(method)
    type <- match.arg(type, several.ok = TRUE)
    method_data <- read_method_data(
        simulation_id = simulation_id,
        root = root,
        method = method
    )

    out <- lapply(type, .save_method_curves, simulation_id = simulation_id,
                  root = root, method = method, method_data = method_data,
                  x0 = x0)
    names(out) <- type

    invisible(out)
}

.save_method_curves <- function(type, simulation_id, root, method,
                                method_data, x0) {
    curve_fn <- switch(type,
                       roc = make_average_roc_curve,
                       fdr = make_average_fdr_curve)
    curve <- curve_fn(method_data = method_data, x0 = x0)
    path <- get_curve_path(simulation_id, root, method, type = type)
    message(glue::glue("Writing {toupper(type)} curve to {path}..."))
    write.table(curve, path)
    curve
}

#' Read in the ROC or FDR curve data for a given method
#' @inheritParams read_method_data
#' @export
read_method_curves <- function(
        simulation_id,
        root = "out/simulation_studies",
        method = c("edgeR",
                   "DESeq2",
                   "limma",
                   "ngstan"),
        type = c("roc", "fdr")
) {
    method <- match.arg(method, several.ok = TRUE)
    names(method) <- method
    type <- match.arg(type, several.ok = TRUE)
    names(type) <- type
    lapply(type, .read_method_curves, simulation_id = simulation_id,
           root = root, method = method)
}

.read_method_curves <- function(method, simulation_id, root, type) {
    if (length(method) > 1) {
        curves <- lapply(method, .read_method_curves,
                         simulation_id = simulation_id,
                         root = root, type = type)
        return(dplyr::bind_rows(curves, .id = "method"))
    }
    read.table(
        get_curve_path(simulation_id, root, method, type)
    )
}

#' Create plots of ROC and FDR curves
#' @inheritParams read_method_data
#' @export
plot_curves <- function(
        simulation_id,
        root = "out/simulation_studies",
        method = c("edgeR",
                   "DESeq2",
                   "limma",
                   "ngstan"),
        type = c("roc", "fdr")
) {
    method <- match.arg(method, several.ok = TRUE)
    type <- match.arg(type, several.ok = TRUE)
    curve_data <- read_method_curves(simulation_id = simulation_id,
                                     root = root, method = method, type = type)
    purrr::imap(curve_data, ~ {.plot_curves(.x, .y)})
}

.plot_curves <- function(x, type = c("roc", "fdr")) {
    xaes <- ggplot2::sym(switch(type,
                                roc = "fpr",
                                fdr = "fdr"))
    yaes <- ggplot2::sym(switch(type,
                                roc = "average_tpr",
                                fdr = "tfdr"))
    ggplot2::ggplot(x, ggplot2::aes(!!xaes, !!yaes, color = method)
    ) + ggplot2::geom_line()
}

#' Create and save ROC and FDR curves for a simulation study
#' @inheritParams read_method_data
#' @export
save_curve_plots <- function(
        simulation_id,
        root = "out/simulation_studies",
        method = c("edgeR",
                   "DESeq2",
                   "limma",
                   "ngstan"),
        type = c("roc", "fdr")
) {
    method <- match.arg(method, several.ok = TRUE)
    type <- match.arg(type, several.ok = TRUE)
    plots <- plot_curves(simulation_id = simulation_id, root = root,
                         method = method, type = type)
    purrr::imap(plots, .save_curve_plots, simulation_id = simulation_id, root = root)
}

.save_curve_plots <- function(plot, plot_name, simulation_id, root) {
    fn <- file.path(root, simulation_id, glue::glue("{simulation_id}_{plot_name}.pdf"))
    message(glue::glue("Saving {plot_name} plot to {fn}..."))
    ggplot2::ggsave(filename = fn, plot = plot)
}
