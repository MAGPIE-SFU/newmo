library(dplyr)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)


parameter_generator <- function(
    how_numerical = "SSE",
    optimizer = "hjkb",
    vaccine_type1 = "PCV7",
    vaccine_type2 = "PCV13",
    use_switch = FALSE,
    how_gene_order = "default"
) {

    if (how_numerical == "CC_MAGIC") {
        # Read optimizer output
        parameters <- readRDS(common$CC_MAGIC)
    } else {
        # Read optimizer output
        parameters <- readRDS(common$FITTED_PARAMETERS_PATH)

        # Select parameter set by loss function and optimizer
        idx <- which(
            parameters$loss == how_numerical &
            parameters$method == optimizer
        )

        if (length(idx) != 1) {
            stop(
                "Expected exactly one parameter set for loss = '",
                how_numerical, "' and optimizer = '", optimizer,
                "'; found ", length(idx), "."
            )
        }
    }

    modelpars <- as.list(
        parameters[idx, setdiff(names(parameters), "loss"), drop = FALSE]
    )

    # ------------------------------------------------
    # Vaccine related input
    # ------------------------------------------------
    modelpars$tswitch <- if (use_switch) common$T_SWITCH else 0

    # Each argument may be either a name looked up in common$VACCINE_TYPES
    # or a raw character vector of strain labels.
    resolve_vaccine <- function(spec, arg_name) {
        if (is.character(spec) && length(spec) == 1L && !is.na(spec) &&
            spec %in% names(common$VACCINE_TYPES)) {
            return(common$VACCINE_TYPES[[spec]])
        }
        if (is.character(spec) && length(spec) > 0L && !anyNA(spec)) {
            return(spec)
        }
        stop(
            arg_name, " must name a vaccine in common$VACCINE_TYPES or be ",
            "a non-empty character vector of strain labels"
        )
    }

    modelpars$vaccine1 <- resolve_vaccine(vaccine_type1, "vaccine_type1")
    modelpars$vaccine2 <- resolve_vaccine(vaccine_type2, "vaccine_type2")

    # ------------------------------------------------
    # Gene-order parameterization
    # ------------------------------------------------
    locus_by_year <- read_csv(common$LOCUS_FILE_PATH, show_col_types = FALSE)
    locus_by_year <- filter(locus_by_year, in_model)
    # note we rely on the order of the loci
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


    if (identical(how_gene_order, "default")) {
        modelpars$gene_order <- NULL
    } else if (
        is.character(how_gene_order) &&
        length(how_gene_order) == 1
    ) {
        modelpars$gene_order <-  switch(
            how_gene_order,
            delta      = order_by("delta"),
            z          = order_by("z"), # used to be z; number of sds that the freq mvoed, under a simple model for the sd
            lograt = order_by("lograt"),
            logitdiff = order_by("logitdiff"),
            varlogitpre = order_by("varlogitpre"),
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
