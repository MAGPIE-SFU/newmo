#' Serotypes covered by each pneumococcal vaccine
#'
#' A named list of character vectors, one per vaccine.
#' Labels must match the serotype labels used in the model data
#' (check "15B"/"15C" and "20"/"20A").
#' Renamed 15B in PCV21 and PCV23 to 15B/C due to other parts of our model
#'
#' @export
VACCINE_TYPES <- list(
  PCV7 = c(
    "4", "6B", "9V", "14", "18C", "19F", "23F"
  ),
  PCV13 = c(
    "1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C",
    "19A", "19F", "23F"
  ),
  PCV15 = c(
    "1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C",
    "19A", "19F", "22F", "23F", "33F"
  ),
  PCV21 = c(
    "3", "6A", "7F", "8", "9N", "10A", "11A", "12F", "15A", "15B/C",
    "16F", "17F", "19A", "20A", "22F", "23A", "23B", "24F", "31", "33F",
    "35B"
  ),
  PCV23 = c(
    "1", "2", "3", "4", "5", "6B", "7F", "8", "9N", "9V",
    "10A", "11A", "12F", "14", "15B/C", "17F", "18C", "19A", "19F", "20",
    "22F", "23F", "33F"
  )
)

#' Get the serotypes of one vaccine
#'
#' @param name A vaccine name, e.g. "PCV13".
#' @return A character vector of serotype labels.
#' @export
vaccine_serotypes <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single vaccine name", call. = FALSE)
  }
  if (!name %in% names(VACCINE_TYPES)) {
    stop("Unknown vaccine '", name, "'. Available: ",
         paste(names(VACCINE_TYPES), collapse = ", "), call. = FALSE)
  }
  VACCINE_TYPES[[name]]
}
