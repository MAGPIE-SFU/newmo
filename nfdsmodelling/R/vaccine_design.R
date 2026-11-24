# ============================================================
#  vaccine_design.R
#  NFDS-based vaccine design (was explore-vaccine-fns.R and the
#  vaccine design script): score a vaccine by the IPD it leaves,
#  and invent vaccines by IPD fraction, invasiveness, a greedy
#  approach, and an exhaustive search.
# ============================================================

#' Vaccine as a single string
#' @param vaccine Character vector of serotypes.
#' @export
compress_vstring <- function(vaccine) { paste(vaccine, collapse="-")}

#' Predicted IPD of a simulation
#'
#' @param simulation A simulation (rows = groups, columns year.1 ...).
#' @param inv Invasiveness lookup (serotype, v_child).
#' @param group_lookup Group definitions (group, serotype).
#' @param years Simulation years to average over, e.g. 6:10.
#' @return A list with `total` (summed IPD) and `seroIPD` (IPD per serotype).
#' @importFrom dplyr %>% select inner_join mutate
#' @export
getIPD <- function(simulation, inv, group_lookup, years = 6) {
  cols <- paste0("year.", years)
  missing <- setdiff(cols, names(simulation))
  if (length(missing)) stop("Years not in simulation: ",
                            paste(missing, collapse = ", "))

  # simulation needs the serotype, right now the row names are the groups
  simulation$serotype = group_lookup$serotype[match( rownames(simulation), group_lookup$group)]
  # mean IPD in the years we are using
  simulation$xmean <- rowMeans(simulation[, cols, drop = FALSE])

  dat = simulation %>% select(serotype, xmean) %>% inner_join( inv, by = "serotype") %>%
    mutate(IPD = xmean*v_child)

  totalIPD = sum(dat$IPD, na.rm = T)
  return(list(total = totalIPD, seroIPD = dat %>% select(serotype, IPD)))
}

#' Score a vaccine by the IPD it leaves
#'
#' Simulates with `vseros` as the only vaccine (vaccine1 = vaccine2) and
#' returns the total IPD, averaged over `years`.
#'
#' @param modelpars Model parameters for [simulate_nfds()].
#' @param vseros Serotypes in the vaccine.
#' @param samodel The model structure.
#' @inheritParams getIPD
#' @export
score_vaccine <- function(modelpars, vseros, samodel, inv, group_lookup, years = 6) {
  modelpars$vaccine1=vseros
  modelpars$vaccine2 = vseros
  mysim = simulate_nfds(samodel, modelpars, nyears=max(years))
  ipd  <- getIPD(mysim, inv, group_lookup, years = years)$total # only cumulative score
  ipd
}

#' Pre-vaccine serotype composition, invasiveness and vaccine membership
#'
#' Serotype-level carriage fraction in the inferred pre-vaccine population,
#' child invasiveness, their product (the IPD contribution), and the first
#' vaccine in `vaccines` that contains each serotype ("NVT" if none).
#'
#' @param xt Grouped x(t), with a serotype column.
#' @param inv Invasiveness lookup (serotype, v_child).
#' @param vaccines Named list of vaccines, in priority order (e.g. PCV7, PCV13,
#'   PCV15, PCV20, PCV23): a serotype gets the first vaccine it is in.
#' @param prevax_years Years for the pre-vaccine composition.
#' @return One row per serotype: serotype, frac, v_child, ipdfrac, ipdshare, vax.
#' @importFrom dplyr %>% select filter group_by summarise inner_join mutate
#' @export
vaccine_candidate_data <- function(xt, inv, vaccines, prevax_years = 2000:2009) {

  sero_inv <- inv %>%
    select(serotype, v_child) %>%
    filter(!is.na(v_child))

  #data on the composition of the inferred carriage population before vaccination
  compodat <- xt %>%
    filter(year %in% prevax_years) %>%        # e.g. pre-vaccine only
    group_by(year, serotype) %>%
    summarise(frac = sum(frac), .groups = "drop") %>%   # over groups
    group_by(serotype) %>%
    summarise(frac = mean(frac), .groups = "drop")  %>%
    inner_join(sero_inv, by = "serotype") %>%
    mutate(ipdfrac = v_child*frac)%>% # does not sum to 1, but comparable from simulation to simulation
    mutate(ipdshare = ipdfrac / sum(ipdfrac)) # normalized, sums to 1

  # the first vaccine (in the given order) that contains each serotype
  vax_of <- function(s) {
    for (v in names(vaccines)) if (s %in% vaccines[[v]]) return(v)
    "NVT"
  }
  compodat %>% mutate(vax = vapply(serotype, vax_of, character(1)))
}

#' Invent vaccines
#'
#' Method 1 ("ipdfrac", "mostinvasive"): the top `direct_size` serotypes by
#' their contribution to IPD in the inferred pre-vaccine population, or by
#' child invasiveness. Method 2 ("greedy"): starting from `base_vaccine`, add
#' the candidate that lowers IPD most, one at a time, up to `target_size`.
#' Method 3 ("exhaustive"): every combination of candidates added to
#' `base_vaccine` up to `target_size`.
#'
#' Candidates for methods 2 and 3: every serotype in one of
#' `candidate_vaccines`, plus the NVT serotypes in the top `n_extra` by IPD
#' contribution or by child invasiveness, minus `base_vaccine`.
#'
#' @param modelpars Model parameters for [simulate_nfds()] (vaccines are set per score).
#' @param samodel The model structure.
#' @param xt Grouped x(t), with a serotype column.
#' @param inv Invasiveness lookup.
#' @param group_lookup Group definitions.
#' @param base_vaccine Serotypes every greedy / exhaustive vaccine starts from (PCV7).
#' @param candidate_vaccines Named list of vaccines whose serotypes are candidates,
#'   in priority order (see [vaccine_candidate_data()]).
#' @param prevax_years Years for the inferred pre-vaccine composition.
#' @param direct_size Valency of the method 1 vaccines.
#' @param target_size Valency of the greedy and exhaustive vaccines.
#' @param n_extra Number of top NVT serotypes (by IPD contribution, and by
#'   invasiveness) added to the candidates.
#' @param years Simulation years to average the IPD over.
#' @param n_cores Cores for the exhaustive search (> 1 uses mclapply, Linux/macOS only).
#' @return One row per vaccine: ipd, vaccine (string), method.
#' @importFrom dplyr %>% select filter slice_max pull tibble add_row arrange bind_rows
#' @importFrom utils combn
#' @export
design_vaccines <- function(modelpars, samodel, xt, inv, group_lookup,
                            base_vaccine, candidate_vaccines,
                            prevax_years = 2000:2009,
                            direct_size = 13, target_size = 15, n_extra = 5,
                            years = 6:10, n_cores = 1) {

  compodat <- vaccine_candidate_data(xt, inv, candidate_vaccines, prevax_years)

  # ---- method 1: use just the estimated contribution to the IPD population in our inferred dataset

  # the top contributors to IPD in our inferred pre-vaccine population
  top13 <- compodat %>%
    slice_max(ipdfrac, n = direct_size) %>%
    pull(serotype) %>%
    as.character() # as it happens this does include PCV7

  vax_pre_ipdfrac = top13 # the vaccine you would choose in method 1.

  # what about the top most invasive strains? (not fraction contribs)?
  top13car <- compodat %>%
    slice_max(v_child, n = direct_size) %>%
    pull(serotype) %>%
    as.character()

  direct = tibble(vaccine = compress_vstring(vax_pre_ipdfrac),
                  ipd = score_vaccine(modelpars, vax_pre_ipdfrac, samodel, inv, group_lookup, years=years),
                  method = "ipdfrac")
  direct = rbind(direct,
                 tibble(vaccine = compress_vstring(top13car),
                        ipd = score_vaccine(modelpars, top13car, samodel, inv, group_lookup, years=years),
                        method = "mostinvasive"))

  # ---- candidate set: serotypes in the candidate vaccines, plus the top NVT
  # serotypes by IPD contribution or by invasiveness
  in_vax <- compodat %>% filter(vax != "NVT") %>% pull(serotype)

  extra <- compodat %>%
    filter(vax == "NVT") %>%
    filter(rank(-ipdfrac) <= n_extra | rank(-v_child) <= n_extra) %>%
    pull(serotype)

  all_candidates = union(in_vax, extra)
  candidates <- setdiff(all_candidates, base_vaccine)

  n_add <- target_size - length(base_vaccine)
  message(length(candidates), " candidates; the exhaustive search scores ",
          format(choose(length(candidates), n_add), big.mark = ","), " vaccines")

  # ---- method 2: greedy .
  # explanation: of the candidate serotypes (according to their estimated invasiveness in
  # children or their IPD contribution or their existence in other PCV vaccines aimed at
  # children), choose the best one to add to PCV7. Then choose the best one to add to
  # c(PCV, firstone). Etc until we get to a target_size-valent vaccine.

  current   <- base_vaccine
  remaining <- candidates
  trace     <- tibble(step = integer(), vaccine = character(),
                      ipd = double())

  tmp = base_vaccine
  while (length(current) < target_size && length(remaining) > 0) {

    scores <- vapply(remaining,
                     function(s) score_vaccine(modelpars, c(current, s), samodel, inv, group_lookup, years = years),
                     FUN.VALUE = 1)

    best <- remaining[which.min(scores)]

    trace <- add_row(trace,
                     step     = length(current) - length(base_vaccine) + 1,
                     vaccine    = compress_vstring(c(tmp,best)),
                     ipd      = min(scores))
    tmp = c(tmp, best)
    current   <- c(current, best)
    remaining <- setdiff(remaining, best)

    message(sprintf("step %d: added %s, IPD %.1f",
                    nrow(trace), best, min(scores)))
  }

  # ---- method 3: exhaustive search
  # all possible (n candidates choose n_add) combinations of the candidates
  combos <- combn(candidates, n_add, simplify = FALSE)   # combos[[1]] is a character vector of things added to PCV7

  # for each one, add those to PCV7 and compute the score
  score_combo <- function(s) score_vaccine(modelpars, c(base_vaccine, s), samodel, inv, group_lookup, years = years)
  scores <- if (n_cores > 1) {
    unlist(parallel::mclapply(combos, score_combo, mc.cores = n_cores))
  } else {
    vapply(combos, score_combo, FUN.VALUE = 1)
  }

  # collect the results
  exhaustive_results <- tibble(
    id    = seq_along(combos),
    vaccine = vapply(combos,
                     function(s) paste(c(base_vaccine, s), collapse = ","),
                     FUN.VALUE = character(1)),
    ipd   = scores
  ) %>%
    arrange(ipd)

  # join all the results together
  trace$method= "greedy"
  exhaustive_results$method = "exhaustive"

  bind_rows(
    trace               %>% select(ipd, vaccine, method),
    direct              %>% select(ipd, vaccine, method),
    exhaustive_results  %>% select(ipd, vaccine, method)
  )
}

#' Invasiveness and IPD contribution of each serotype, by vaccine
#'
#' Illustrates which serotypes are in which vaccines, and their (1) child
#' invasiveness and (2) contribution to pre-vaccine IPD. Shows why the
#' candidates for the vaccine design are chosen.
#'
#' @param compodat From [vaccine_candidate_data()].
#' @param vax_cols Named colours per vaccine (and "NVT").
#' @return A patchwork of two plots.
#' @importFrom ggplot2 ggplot aes geom_segment scale_colour_manual scale_x_log10
#'   expansion coord_cartesian labs theme_bw
#' @importFrom patchwork plot_layout
#' @export
plot_vaccine_candidates <- function(compodat,
                                    vax_cols = c(
                                      PCV7  = "#440154",   # dark purple
                                      PCV13 = "#414487",
                                      PCV15 = "#2a788e",
                                      PCV20 = "#22a884",
                                      PCV23 = "#7ad151",   # bright green
                                      NVT   = "grey70"
                                    )) {

  V_FLOOR   <- 10^floor(log10(min(compodat$v_child)))
  # bars on a log axis can't start at 0, so they start at the next power of 10 below the smallest value
  IPD_FLOOR <- 10^floor(log10(min(compodat$ipdfrac[compodat$ipdfrac > 0])))

  # invasiveness, for reference
  p_inv <- ggplot(compodat, aes(y = serotype, colour = vax)) +
    geom_segment(aes(x = V_FLOOR, xend = v_child, yend = serotype),
                 linewidth = 4, lineend = "butt") +
    scale_colour_manual(values = vax_cols, name = NULL) +
    scale_x_log10(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(V_FLOOR, max(compodat$v_child))) +
    labs(x = "invasiveness (children)", y = NULL) +
    theme_bw()

  p_ipdfrac <- ggplot(compodat, aes( y = serotype, colour = vax)) +
    geom_segment(aes(x = IPD_FLOOR, xend = ipdfrac, yend = serotype),
                 linewidth = 4, lineend = "butt") +
    scale_x_log10(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(IPD_FLOOR, max(compodat$ipdfrac))) +
    scale_colour_manual(values = vax_cols, name = NULL) +
    labs(x = "fraction of inferred pre-vaccine ipd", y = NULL) +
    theme_bw()

  p_inv + p_ipdfrac + plot_layout(guides = "collect")
}