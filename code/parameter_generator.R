parameter_generator <- function(
    how_numerical = "CC_MAGIC",
    vaccine_type = "PCV13",
    how_gene_order = "default"
) {
    # Read optimizer output
    parameters <- readRDS("data/SouthAfrica/fitted_parameters.rds")

    # Select parameter set
    idx <- match(how_numerical, parameters$loss)

    if (is.na(idx)) {
        stop("Unknown parameter set: ", how_numerical)
    }

    modelpars <- as.list(
        parameters[idx, setdiff(names(parameters), "loss"), drop = FALSE]
    )

    # Vaccine type
    if (identical(vaccine_type, "PCV13")) {
        modelpars$vaccine <- c(
            "1", "3", "4", "5", "6A", "6B", "7F",
            "9V", "14", "18C", "19A", "19F", "23F"
        )
    } else if (identical(vaccine_type, "PCV7")) {
        modelpars$vaccine <- c(
            "4", "6B", "9V", "14", "18C", "19F", "23F"
        )
    } else if (
        is.character(vaccine_type) &&
        length(vaccine_type) > 0 &&
        !anyNA(vaccine_type)
    ) {
        modelpars$vaccine <- vaccine_type
    } else {
        stop(
            "vaccine_type must be 'PCV13', 'PCV7', ",
            "or a non-empty character vector"
        )
    }

    # Gene-order parameterization
    if (identical(how_gene_order, "default")) {
        modelpars$gene_order <- NULL
    } else if (
        is.character(how_gene_order) &&
        length(how_gene_order) == 1
    ) {
        switch(
            how_gene_order,
            logit_diff = stop("logit_diff is not implemented yet"),
            delta = stop("delta is not implemented yet"),
            zett = stop("zett is not implemented yet"),
            log_ratio = stop("log_ratio is not implemented yet"),
            stop("Unknown gene-order method: ", how_gene_order)
        )
    } else if (
        is.numeric(how_gene_order) ||
        is.list(how_gene_order)
    ) {
        modelpars$gene_order <- how_gene_order
    } else {
        stop("Invalid how_gene_order")
    }

    modelpars
}
