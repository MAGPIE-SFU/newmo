common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)


parameter_generator <- function(
    how_numerical = "CC_MAGIC",
    vaccine_type = "PCV13",
    how_gene_order = "default"
) {
    # Read optimizer output
    parameters <- readRDS(common$FITTED_PARAMETERS_PATH)

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
        modelpars$vaccine <- common$VACCINE_TYPES[["PCV13"]]
    } else if (identical(vaccine_type, "PCV7")) {
        modelpars$vaccine <- common$VACCINE_TYPES[["PCV7"]]
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

    locus_by_year <- read_csv(common$LOCUS_FILE_PATH, show_col_types = FALSE)
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
            z     = order_by("z"), # used to be z; number of sds that the freq mvoed, under a simple model for the sd
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
