parameter_generator <- function(
    how_numerical = "CC_MAGIC",
    vaccine_type = "PCV13",
    how_gene_order = "default",
    locus_file = "data/SouthAfrica/processed-inputs/locus_freqs_deltas.csv"
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
    
    locus_by_year <- read_csv(locus_file, show_col_types = FALSE)
    locus_by_year <- filter(locus_by_year, in_model) # note we rely on the order of the loci
    # being the same as in the model's G matrix, or we load the samodel here. 
    # this is checked in the file that creates the SA nfds model (that the locus order is correct) 

    order_by <- function(col, decreasing = FALSE) {
        if (is.null(locus_by_year))
            stop("need locus_by_year information")
        if (!col %in% names(locus_by_year))
            stop("locus_by_year has no column '", col, "'")
        v <- locus_by_year[[col]]
        order(v, decreasing = decreasing)
    }
    
    # Gene-order parameterization
    if (identical(how_gene_order, "default")) {
        modelpars$gene_order <- NULL
    } else if (
        is.character(how_gene_order) &&
        length(how_gene_order) == 1
    ) {
        modelpars$gene_order <-  switch(
            how_gene_order,
            delta      = order_by("delta"),
            z          = order_by("z"),
            lograt  = order_by("lograt"),
            logitdiff = order_by("logitdiff"),
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
