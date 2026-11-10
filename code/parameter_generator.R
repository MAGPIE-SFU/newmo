# ============================================================
#  parameter_generator.R
#  Build the modelpars list used by simulate() from saved fits.
# ============================================================

library(dplyr)
library(readr)

# tmp loadign the local package to get some things working...
pkgload::load_all(here::here("nfdsmodelling"), quiet = TRUE)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

# Columns that identify a fit rather than parameterize the model
FIT_ID_COLS <- c("fit_id", "loss", "method", "start_id")

# ------------------------------------------------------------
# Vaccines
# ------------------------------------------------------------

# A spec is either a name in common$VACCINE_TYPES or a raw character vector of
# strain labels. An unknown single string is an ERROR (not silently treated as
# a strain label that matches nothing, which would mean no vaccination).
resolve_vaccine <- function(spec, arg_name) {
  ok_type <- is.character(spec) && length(spec) > 0L && !anyNA(spec)
  if (!ok_type) {
    stop(arg_name, " must be a vaccine name or a non-empty character vector of strain labels")
  }

  if (length(spec) == 1L) {
    return(nfdsmodelling::VACCINE_TYPES[[spec]])
  }

  spec
}

# ------------------------------------------------------------
# Loading fitted parameters
# ------------------------------------------------------------

load_fitted_row <- function(how_numerical, optimizer, use_switch, fit_id, path) {
  if (is.null(path)) {
    path <- if (use_switch) common$FITTED_PARAMETERS_PATH_SWITCH
            else            common$FITTED_PARAMETERS_PATH
  }
  parameters <- readRDS(path)

  if (!is.null(fit_id)) {
    if (!"fit_id" %in% names(parameters)) {
      stop("'", basename(path), "' has no fit_id column")
    }
    idx <- which(parameters$fit_id == fit_id)
  } else {
    idx <- which(parameters$loss == how_numerical & parameters$method == optimizer)
  }

  if (length(idx) != 1) {
    stop(
      "Expected exactly one parameter set in '", basename(path), "' for ",
      if (!is.null(fit_id)) paste0("fit_id = '", fit_id, "'")
      else paste0("loss = '", how_numerical, "' and optimizer = '", optimizer, "'"),
      "; found ", length(idx), ". ",
      "If the file holds several starts, use a best-per-fit file or pass fit_id."
    )
  }

  parameters[idx, , drop = FALSE]
}

# ------------------------------------------------------------
# Gene order
# ------------------------------------------------------------

GENE_ORDER_COLUMNS <- c("delta", "z", "lograt", "logitdiff", "varlogitpre")

load_locus_by_year <- function() {
  locus_by_year <- read_csv(common$LOCUS_FILE_PATH, show_col_types = FALSE)
  # NB: relies on the locus order matching the model's G matrix. That is
  # checked in the script that creates the SA nfds model.
  filter(locus_by_year, in_model)
}

validate_gene_order <- function(x, n_loci) {
  x <- unlist(x)
  if (!is.numeric(x) || length(x) != n_loci || !setequal(x, seq_len(n_loci))) {
    stop("gene_order must be a permutation of 1..", n_loci)
  }
  x
}

resolve_gene_order <- function(how_gene_order, decreasing = FALSE) {
  if (identical(how_gene_order, "default")) return(NULL)

  locus_by_year <- load_locus_by_year()
  n_loci <- nrow(locus_by_year)

  if (is.character(how_gene_order) && length(how_gene_order) == 1L) {
    if (!how_gene_order %in% GENE_ORDER_COLUMNS) {
      stop("Unknown gene-order method: ", how_gene_order,
           " (use 'default', one of ", paste(GENE_ORDER_COLUMNS, collapse = ", "),
           ", or a permutation)")
    }
    if (!how_gene_order %in% names(locus_by_year)) {
      stop("locus_by_year has no column '", how_gene_order, "'")
    }
    return(order(locus_by_year[[how_gene_order]], decreasing = decreasing))
  }

  if (is.numeric(how_gene_order) || is.list(how_gene_order)) {
    return(validate_gene_order(how_gene_order, n_loci))
  }

  stop("Invalid how_gene_order")
}

# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------

# NOTE: the vaccine defaults must match FIXED_PARS in the fitting functions,
# otherwise parameters fitted for one vaccine setup are reused with another.
parameter_generator <- function(
    how_numerical = "SSE",
    optimizer = "hjkb",
    vaccine_type1 = "PCV7",
    vaccine_type2 = "PCV13",
    use_switch = FALSE,
    how_gene_order = "default",
    gene_order_decreasing = FALSE,
    fit_id = NULL,     # pick one specific fit, e.g. "weighted_SSE__nmkb__round_05"
    path = NULL        # override the parameter file
) {

  if (how_numerical == "CC_MAGIC") {
    # Baseline = the default parameters, built from common.R (no file needed)
    modelpars <- common$DEFAULT_PARS
    modelpars$wlow  <- modelpars$whigh * modelpars$wfrac
    modelpars$wfrac <- NULL
    modelpars$v2    <- modelpars$v1
  } else {
    row <- load_fitted_row(how_numerical, optimizer, use_switch, fit_id, path)

    # The file should have been fitted with the same switch setting
    if ("tswitch" %in% names(row) && (row$tswitch > 0) != use_switch) {
      warning("use_switch = ", use_switch, " but the loaded parameters were fitted ",
              if (row$tswitch > 0) "with" else "without", " a vaccine switch.")
    }

    modelpars <- as.list(row[, setdiff(names(row), FIT_ID_COLS), drop = FALSE])
  }

  # Vaccines
  modelpars$tswitch  <- if (use_switch) common$T_SWITCH else 0
  modelpars$vaccine1 <- resolve_vaccine(vaccine_type1, "vaccine_type1")
  modelpars$vaccine2 <- resolve_vaccine(vaccine_type2, "vaccine_type2")

  # Gene order (NULL means the default order inside simulate())
  modelpars$gene_order <- resolve_gene_order(how_gene_order, gene_order_decreasing)

  modelpars
}


load_best_fit_pars <- function(path = common$OVERALL_BEST_FITS_PATH,
                               use_switch = TRUE, ...) {
  best <- readRDS(path)
  if (nrow(best) != 1) {
    stop("Expected one row in '", basename(path), "', found ", nrow(best))
  }

  parameter_generator(
    how_numerical = best$loss,
    optimizer     = best$method,
    fit_id        = best$fit_id,
    path          = path,
    use_switch    = use_switch,
    ...
  )
}

