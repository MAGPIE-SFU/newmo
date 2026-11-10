# this pulls the reworked data together and creates the
# NFDS model structure (a list)

#' Build the NFDS model structure
#'
#' @param G Grouped G for the NFDS model, e.g. `build_gnfds()$G`.
#' @param xt Grouped x(t), e.g. from [reweight_carriage_uncertainty()].
#' @param locus_by_year Locus statistics, e.g. `build_gnfds()$locus_by_year`.
#' @param x0_years Years averaged for the starting composition x0.
#' @return A list with G, st (serotype of each row of G), x0 and el
#'   (equilibrium locus frequencies).
#' @importFrom dplyr %>% filter group_by summarise
#' @importFrom stats setNames
#' @export
build_samodel <- function(G, xt, locus_by_year, x0_years) {

  x0 <- xt %>%
    filter(year %in% x0_years) %>%
    group_by(group) %>%
    summarise(frac = mean(frac), .groups = "drop")

  x0 <- setNames(x0$frac, x0$group)[rownames(G)]   # align to G's rows
  stopifnot(!anyNA(x0), abs(sum(x0) - 1) < 1e-8)

  xt$serotype <- ifelse(xt$group == "rare_NVT", "rare_NVT",
                        sub("_[^_]+$", "", xt$group))

  st_year <- min(xt$year)
  stopifnot(identical(rownames(G), xt$group[xt$year == st_year]))
  stopifnot(identical(colnames(G), locus_by_year$locus[locus_by_year$in_model]))
  # this matters in the parameter generator now: otherwise that file needs to load the model
  #  if you want use loci selected on another
  # criterion (z, or logratio, or delta, or whatever) then you can see that in the
  # parameter generator now

  list(G  = G,
       st = xt$serotype[xt$year == st_year], # serotype for each strain (ie each row of G)
       x0 = x0,  # use mean pre-vaccine and just fit to post-vaccine?
       el = as.numeric(locus_by_year$pre2009[locus_by_year$in_model]) # equilibrium frequencies
  )
}
