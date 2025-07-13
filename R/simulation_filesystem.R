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
        docfile <- file.path(root, "study_index.csv")
        docexists <- file.exists(docfile)
        if (docexists) append <- TRUE
        else append <- FALSE
        desc <- readline("Enter a one-line description: ")
        row <- data.frame(simulation_id = simulation_id, description = desc)
        message(glue::glue("Writing description to {docfile}"))
        write.csv(row, file = docfile, append = append, row.names = FALSE)
    }
    return(path)
}
