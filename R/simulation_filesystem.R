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
make_sim_dir <- function(simulation_id, root = "out/simulation_studies") {
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
        if (docexists) append <- TRUE
        desc <- readline("Enter a one-line description: ")
        row <- data.frame(simulation_id = simulation_id, description = desc)
        message(glue::glue("Writing description to {docfile}"))
        write.table(row, file = docfile, append = append, row.names = FALSE)
    }
    return(path)
}

#' Set up a directory to contain all the files for a single simulated dataset
#' within a simulation study
#' @param dataset_id The dataset identifier
#' @inheritParams make_sim_dir
#' @export
make_dataset_dir <- function(dataset_id, simulation_id, root = "out/simulation_studies") {
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
    if (docexists) append <- TRUE
    row <- data.frame(dataset_id = dataset_id)
    write.table(row, file = docfile, append = append, row.names = FALSE)
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
