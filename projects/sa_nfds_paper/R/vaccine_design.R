# Model parameters from a fit row, for scoring vaccines (vaccines are set per score)
fit_modelpars <- function(fit_row, tswitch) {
  modelpars <- as.list(fit_row[c("v1", "v2", "sigma", "portion", "whigh", "wlow", "migr")])
  modelpars$tswitch  <- tswitch
  modelpars$vaccine1 <- VACCINE_TYPES[["PCV7"]]
  modelpars$vaccine2 <- VACCINE_TYPES[["PCV13"]]
  modelpars
}
