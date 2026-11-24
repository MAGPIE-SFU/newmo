# ============================================================
# Simulates an expected IPD population from a simulated carriage pop
# of n_carriage per year, includes uncertainty propagation from invasiveness
# estimates. plot_projected_ipd_check() plots it against the observed IPD.
# ============================================================

#' Simulate projected IPD with uncertainty from invasiveness estimates
#'
#' Scales a simulated carriage composition to `n_carriage` per year, splits
#' each group over its serotypes, and multiplies by child invasiveness drawn
#' `n_reps` times from its log-normal distribution. Simulated carriage is held
#' fixed. `ipd_expected` is a Poisson mean, not a realised case count.
#'
#' @param sim A simulation from [run_model()] / [simulate_nfds()]: rows =
#'   groups, columns year.1 ... year.n.
#' @param samodel The model structure the simulation was run with.
#' @param inv Invasiveness lookup, e.g. from [build_invasiveness_lookup()].
#' @param group_lookup Group definitions, e.g. `reweight_carriage()$group_lookup`.
#' @param meta_clean Cleaned metadata, e.g. from [clean_metadata()].
#' @param non_typeable Serotype labels with zero invasiveness.
#' @param n_reps Number of Monte Carlo draws.
#' @param n_carriage Simulated carriage population per year.
#' @param seed Random seed, set at the start as in the original script.
#'
#' @return A list with `ipd_summary`, `ipd_yearly`, `ipd_yearly_draws` and
#'   `ipd_serotype_summary`.
#'
#' @importFrom dplyr %>% mutate group_by ungroup summarise transmute if_else
#'   select left_join distinct filter count case_when n first pull arrange
#'   starts_with
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom purrr map_dfr
#' @importFrom tibble rownames_to_column
#' @importFrom stats rlnorm sd var quantile
#' @export
simulate_projected_ipd <- function(sim, samodel, inv, group_lookup, meta_clean,
                                   non_typeable,
                                   n_reps = 1000L,
                                   n_carriage = 100000,
                                   seed = 1) {
  set.seed(seed)

  # 1. Simulated carriage composition ---------------------------------------
  stopifnot(nrow(sim) == length(samodel$st),
            !is.null(rownames(sim)), !anyDuplicated(rownames(sim)))

  sim_long <- sim %>%
    as.data.frame() %>%
    rownames_to_column("group") %>%
    pivot_longer(cols = starts_with("year."),
                 names_to = "year", values_to = "x") %>%
    mutate(group = as.character(group),
           year = as.integer(sub("year.", "", year))) %>%
    group_by(year) %>%
    mutate(frac = x / sum(x), carriage = n_carriage * frac) %>%
    ungroup()
  stopifnot(all(is.finite(sim_long$frac)), all(sim_long$frac >= 0))

  # 2. The SAME child progression-rate sampling scheme as the revised
  #    carriage uncertainty script ------------------------------------------
  inv_dist <- inv %>%
    transmute(serotype = as.character(serotype),
              v_mean = v_child,
              logmean = v_child_logmean,
              logsd = v_child_logsd) %>%
    mutate(v_mean = if_else(serotype %in% non_typeable, 0, v_mean))
  stopifnot(!anyDuplicated(inv_dist$serotype))

  # Each iteration uses one child-rate draw per serotype, held fixed across
  # all simulation years and all groups. The same draw is used for IPD totals,
  # group composition and serotype composition.
  draw_v <- function() {
    z <- inv_dist
    zero_type <- z$serotype %in% non_typeable
    valid <- is.finite(z$logmean) & is.finite(z$logsd) &
      z$logsd >= 0 & is.finite(z$v_mean) & z$v_mean > 0
    invalid <- !zero_type & !valid
    if (any(invalid)) {
      stop("Missing/invalid child log-normal progression-rate parameters for: ",
           paste(z$serotype[invalid], collapse = ", "))
    }
    z$v_sample <- 0
    z$v_sample[valid & !zero_type] <- rlnorm(
      sum(valid & !zero_type),
      meanlog = z$logmean[valid & !zero_type],
      sdlog = z$logsd[valid & !zero_type]
    )
    stopifnot(all(is.finite(z$v_sample)), all(z$v_sample >= 0))
    z %>% select(serotype, v_sample)
  }

  # 3. Map unchanged simulation groups to serotypes ------------------------
  group_lookup <- group_lookup %>%
    mutate(group = as.character(group), serotype = as.character(serotype))
  meta_clean <- meta_clean %>%
    mutate(serotype = as.character(serotype))

  stopifnot(!anyDuplicated(group_lookup[c("serotype", "gpsc")]))

  # Allocate each simulated group among its constituent serotypes.
  # This weighting is only necessary for mixed-serotype groups (e.g. rare_NVT).
  # Use observed child carriage counts across all years, as in the old script.
  serotype_weights <- group_lookup %>%
    distinct(group, serotype) %>%
    left_join(
      meta_clean %>%
        filter(age_group == "child", collection_type == "carriage") %>%
        count(serotype, gpsc, name = "n") %>%
        left_join(group_lookup %>% distinct(serotype, gpsc, group),
                  by = c("serotype", "gpsc"), relationship = "many-to-one") %>%
        group_by(group, serotype) %>%
        summarise(n = sum(n), .groups = "drop"),
      by = c("group", "serotype"), relationship = "one-to-one"
    ) %>%
    mutate(n = replace_na(n, 0)) %>%
    group_by(group) %>%
    mutate(n_serotypes = n(),
           weight = case_when(
             n_serotypes == 1 ~ 1,
             sum(n) > 0 ~ n / sum(n),
             TRUE ~ NA_real_
           )) %>%
    ungroup()

  if (anyNA(serotype_weights$weight)) {
    print(serotype_weights %>% filter(is.na(weight)))
    stop("No observed child-carriage weights for one or more mixed-serotype groups")
  }
  stopifnot(all(abs(serotype_weights %>% group_by(group) %>%
                      summarise(s = sum(weight), .groups = "drop") %>%
                      pull(s) - 1) < 1e-8))

  missing_groups <- setdiff(unique(sim_long$group), unique(serotype_weights$group))
  if (length(missing_groups))
    stop("Unmatched simulated groups: ", paste(missing_groups, collapse = ", "))

  carriage_serotype <- sim_long %>%
    left_join(serotype_weights %>% select(group, serotype, weight),
              by = "group", relationship = "many-to-many") %>%
    mutate(carriage_serotype = carriage * weight)

  missing_v <- setdiff(unique(carriage_serotype$serotype), inv_dist$serotype)
  if (length(missing_v))
    stop("Missing serotypes in invasiveness lookup: ", paste(missing_v, collapse = ", "))

  check_car <- carriage_serotype %>%
    group_by(year) %>%
    summarise(n = sum(carriage_serotype), .groups = "drop")
  stopifnot(all(abs(check_car$n - n_carriage) < 1e-6))

  # 4. Monte Carlo: a SINGLE collection of draws for ALL outputs -----------
  # This propagates uncertainty in estimated rates, holding simulated
  # carriage fixed. `ipd_expected` is a Poisson mean, not a realized case count.
  ipd_serotype_component_draws <- map_dfr(seq_len(n_reps), function(b) {
    v_draw <- draw_v()
    carriage_serotype %>%
      left_join(v_draw, by = "serotype", relationship = "many-to-one") %>%
      mutate(ipd_expected = carriage_serotype * v_sample,
             draw = b) %>%
      select(draw, year, group, serotype, frac, carriage,
             carriage_serotype, ipd_expected)
  })

  # 5. Group-level IPD predictions ------------------------------------------
  ipd_draws <- ipd_serotype_component_draws %>%
    group_by(draw, year, group) %>%
    summarise(frac = first(frac), carriage = first(carriage),
              ipd_expected = sum(ipd_expected), .groups = "drop")

  ipd_summary <- ipd_draws %>%
    group_by(year, group) %>%
    summarise(frac = mean(frac), carriage = mean(carriage),
              ipd_mean = mean(ipd_expected),
              ipd_expected_sd = sd(ipd_expected),
              ipd_sd = sqrt(mean(ipd_expected) + var(ipd_expected)),
              ipd_lower = quantile(ipd_expected, 0.025),
              ipd_upper = quantile(ipd_expected, 0.975),
              .groups = "drop")

  # 6. Total IPD: aggregate INSIDE draws ------------------------------------
  ipd_yearly_draws <- ipd_draws %>%
    group_by(draw, year) %>%
    summarise(total_ipd = sum(ipd_expected), .groups = "drop")

  ipd_yearly <- ipd_yearly_draws %>%
    group_by(year) %>%
    summarise(ipd_mean = mean(total_ipd),
              ipd_expected_sd = sd(total_ipd),
              ipd_sd = sqrt(mean(total_ipd) + var(total_ipd)),
              ipd_lower = quantile(total_ipd, 0.025),
              ipd_upper = quantile(total_ipd, 0.975),
              .groups = "drop")

  # 7. Serotype-level IPD and composition using the SAME draws --------------
  ipd_serotype_draws <- ipd_serotype_component_draws %>%
    group_by(draw, year, serotype) %>%
    summarise(ipd_expected = sum(ipd_expected), .groups = "drop") %>%
    group_by(draw, year) %>%
    mutate(ipd_frac = ipd_expected / sum(ipd_expected)) %>%
    ungroup()

  ipd_serotype_summary <- ipd_serotype_draws %>%
    group_by(year, serotype) %>%
    summarise(frac = mean(ipd_frac), frac_sd = sd(ipd_frac),
              frac_lower = quantile(ipd_frac, 0.025),
              frac_upper = quantile(ipd_frac, 0.975),
              ipd_mean = mean(ipd_expected), .groups = "drop") %>%
    arrange(year, serotype)

  stopifnot(all(abs(ipd_serotype_summary %>%
                      group_by(year) %>% summarise(s = sum(frac), .groups = "drop") %>%
                      pull(s) - 1) < 1e-8))

  # 8. Return outputs (previously written to CSV) ----------------------------
  list(
    ipd_summary          = ipd_summary,
    ipd_yearly           = ipd_yearly,
    ipd_yearly_draws     = ipd_yearly_draws,
    ipd_serotype_summary = ipd_serotype_summary
  )
}

#' Observed vs predicted IPD serotype composition
#'
#' @param ipd_serotype_summary From [simulate_projected_ipd()].
#' @param meta_clean Cleaned metadata (for the observed IPD).
#' @param first_sim_calendar_year Calendar year that corresponds to year.1 in
#'   the simulation.
#' @return A ggplot.
#' @importFrom dplyr %>% filter count group_by mutate ungroup
#' @importFrom ggplot2 ggplot geom_ribbon aes geom_line geom_point facet_wrap labs theme_bw
#' @export
plot_projected_ipd_check <- function(ipd_serotype_summary, meta_clean, first_sim_calendar_year) {

  # Observed vs predicted serotype composition: calendar-year alignment
  observed_ipd <- meta_clean %>%
    filter(collection_type == "IPD") %>%
    count(year, serotype, name = "ipd_count") %>%
    group_by(year) %>%
    mutate(frac = ipd_count / sum(ipd_count)) %>%
    ungroup()

  predicted_ipd <- ipd_serotype_summary %>%
    mutate(year = first_sim_calendar_year + year - 1L, source = "Predicted")

  ggplot() +
    geom_ribbon(data = predicted_ipd,
                aes(x = year, ymin = frac_lower, ymax = frac_upper),
                fill = "#C44E52", alpha = 0.2) +
    geom_line(data = predicted_ipd, aes(year, frac),
              colour = "#C44E52", linewidth = 0.8) +
    geom_point(data = observed_ipd, aes(year, frac),
               colour = "black", size = 1.8) +
    facet_wrap(~serotype, scales = "free_y") +
    labs(x = "Calendar year", y = "Proportion of IPD") +
    theme_bw()
}

#' Evaluate projected IPD under a specified vaccine
#'
#' Simulates pneumococcal carriage following vaccine introduction and
#' projects annual and cumulative serotype-specific invasive pneumococcal
#' disease (IPD), propagating uncertainty in child invasiveness estimates.
#'
#' @param vaccine Character vector of vaccine serotypes.
#' @param modelpars Fitted NFDS model parameters.
#' @param samodel Model structure used for carriage simulation.
#' @param invasiveness_lookup Data frame containing serotype-specific
#'   child invasiveness estimates and log-normal distribution parameters.
#' @param group_lookup Mapping of serotype-GPSC combinations to model groups.
#' @param meta_clean Cleaned pneumococcal metadata.
#' @param non_typeable Character vector of serotypes with zero invasiveness.
#' @param years Number of years following vaccine introduction.
#' @param n_reps Number of Monte Carlo invasiveness draws.
#' @param n_carriage Simulated carriage population per year.
#' @param seed Random seed.
#'
#' @details
#' Mixed-serotype groups are allocated to their constituent serotypes
#' using observed child carriage frequencies.
#'
#' Within each Monte Carlo replicate, one invasiveness value is sampled
#' per serotype and held constant across all simulation years.
#' Cumulative IPD composition is calculated by aggregating expected
#' IPD across years before calculating serotype proportions.
#'
#' Uncertainty intervals reflect uncertainty in invasiveness estimates,
#' conditional on the simulated carriage composition and model parameters.
#'
#' @return A list containing:
#'   \describe{
#'     \item{cumulative}{Tibble of cumulative serotype-specific IPD
#'       burden and composition over the simulation period.}
#'     \item{yearly}{Tibble of annual serotype-specific IPD burden
#'       and composition.}
#'   }
#'
#' @seealso [simulate_nfds()], [score_vaccine()]
#' @export
evaluate_vaccine_ipd <- function(
    vaccine,
    modelpars,
    samodel,
    invasiveness_lookup,
    group_lookup,
    meta_clean,
    non_typeable,
    years = 5L,
    n_reps = 1000L,
    n_carriage = 100000,
    seed = 1L
) {

  set.seed(seed)

  # 1. Simulate carriage under the selected vaccine ----------------------

  modelpars$vaccine1 <- as.character(vaccine)
  modelpars$vaccine2 <- as.character(vaccine)

  sim <- simulate_nfds(
    samodel = samodel,
    modelpars = modelpars,
    nyears = years
  )

  # 2. Annual carriage composition --------------------------------------

  year_cols <- paste0("year.", seq_len(years))

  stopifnot(
    all(year_cols %in% colnames(sim)),
    !is.null(rownames(sim)),
    !anyDuplicated(rownames(sim))
  )

  sim <- as.matrix(sim[, year_cols, drop = FALSE])

  stopifnot(
    all(is.finite(sim)),
    all(sim >= 0),
    all(colSums(sim) > 0)
  )

  # Normalize carriage composition within each year
  sim <- sweep(sim, 2, colSums(sim), "/")

  annual_carriage <- tibble::as_tibble(
    sim * n_carriage,
    rownames = "group"
  ) %>%
    tidyr::pivot_longer(
      cols = -group,
      names_to = "year",
      values_to = "carriage"
    ) %>%
    dplyr::mutate(
      year = as.integer(sub("year.", "", year, fixed = TRUE))
    )

  # 3. Allocate model groups to individual serotypes ---------------------

  group_lookup <- group_lookup %>%
    dplyr::mutate(
      group = as.character(group),
      serotype = as.character(serotype)
    )

  carriage_counts <- meta_clean %>%
    dplyr::filter(
      age_group == "child",
      collection_type == "carriage"
    ) %>%
    dplyr::mutate(
      serotype = as.character(serotype)
    ) %>%
    dplyr::count(serotype, gpsc, name = "n") %>%
    dplyr::left_join(
      group_lookup %>%
        dplyr::distinct(serotype, gpsc, group),
      by = c("serotype", "gpsc"),
      relationship = "many-to-one"
    ) %>%
    dplyr::group_by(group, serotype) %>%
    dplyr::summarise(
      n = sum(n),
      .groups = "drop"
    )

  serotype_weights <- group_lookup %>%
    dplyr::distinct(group, serotype) %>%
    dplyr::left_join(
      carriage_counts,
      by = c("group", "serotype"),
      relationship = "one-to-one"
    ) %>%
    dplyr::mutate(
      n = tidyr::replace_na(n, 0)
    ) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      weight = dplyr::case_when(
        dplyr::n() == 1 ~ 1,
        sum(n) > 0 ~ n / sum(n),
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::ungroup()

  stopifnot(!anyNA(serotype_weights$weight))

  carriage_serotype <- annual_carriage %>%
    dplyr::left_join(
      serotype_weights,
      by = "group",
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      carriage = carriage * weight
    ) %>%
    dplyr::group_by(year, serotype) %>%
    dplyr::summarise(
      carriage = sum(carriage),
      .groups = "drop"
    )

  # Check that annual carriage totals are preserved
  carriage_check <- carriage_serotype %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      total = sum(carriage),
      .groups = "drop"
    )

  stopifnot(
    !anyNA(carriage_serotype$serotype),
    all(abs(carriage_check$total - n_carriage) < 1e-6)
  )

  # 4. Prepare invasiveness distributions -------------------------------

  inv_dist <- invasiveness_lookup %>%
    dplyr::transmute(
      serotype = as.character(serotype),
      v_mean = v_child,
      logmean = v_child_logmean,
      logsd = v_child_logsd
    )

  stopifnot(!anyDuplicated(inv_dist$serotype))

  missing_seros <- setdiff(
    unique(carriage_serotype$serotype),
    inv_dist$serotype
  )

  if (length(missing_seros)) {
    stop(
      "Missing serotypes in invasiveness lookup: ",
      paste(missing_seros, collapse = ", ")
    )
  }

  inv_dist <- inv_dist %>%
    dplyr::filter(
      serotype %in% carriage_serotype$serotype
    ) %>%
    dplyr::arrange(serotype)

  zero_type <- inv_dist$serotype %in% non_typeable

  valid <- is.finite(inv_dist$logmean) &
    is.finite(inv_dist$logsd) &
    inv_dist$logsd >= 0 &
    is.finite(inv_dist$v_mean) &
    inv_dist$v_mean > 0

  if (any(!zero_type & !valid)) {
    stop(
      "Missing or invalid invasiveness parameters for: ",
      paste(
        inv_dist$serotype[!zero_type & !valid],
        collapse = ", "
      )
    )
  }

  # 5. Monte Carlo IPD projections --------------------------------------

  ipd_draws <- purrr::map_dfr(seq_len(n_reps), function(b) {

    # One invasiveness draw per serotype, shared across all years
    v_sample <- numeric(nrow(inv_dist))

    v_sample[!zero_type] <- stats::rlnorm(
      sum(!zero_type),
      meanlog = inv_dist$logmean[!zero_type],
      sdlog = inv_dist$logsd[!zero_type]
    )

    v_draw <- tibble::tibble(
      serotype = inv_dist$serotype,
      v_sample = v_sample
    )

    carriage_serotype %>%
      dplyr::left_join(
        v_draw,
        by = "serotype",
        relationship = "many-to-one"
      ) %>%
      dplyr::transmute(
        draw = b,
        year = year,
        serotype = serotype,
        ipd_expected = carriage * v_sample
      )
  })

  stopifnot(
    all(is.finite(ipd_draws$ipd_expected)),
    all(ipd_draws$ipd_expected >= 0)
  )

  # 6. Annual IPD summary -----------------------------------------------

  ipd_yearly_summary <- ipd_draws %>%
    dplyr::group_by(draw, year) %>%
    dplyr::mutate(
      ipd_frac = ipd_expected / sum(ipd_expected)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(year, serotype) %>%
    dplyr::summarise(
      frac = mean(ipd_frac),
      frac_sd = stats::sd(ipd_frac),
      frac_lower = stats::quantile(ipd_frac, 0.025),
      frac_upper = stats::quantile(ipd_frac, 0.975),

      ipd_mean = mean(ipd_expected),
      ipd_expected_sd = stats::sd(ipd_expected),
      ipd_lower = stats::quantile(ipd_expected, 0.025),
      ipd_upper = stats::quantile(ipd_expected, 0.975),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      vaccine_type = dplyr::if_else(
        serotype %in% vaccine,
        "VT",
        "NVT"
      )
    ) %>%
    dplyr::arrange(year, dplyr::desc(frac))

  # 7. Cumulative IPD summary -------------------------------------------

  ipd_cumulative_summary <- ipd_draws %>%
    dplyr::group_by(draw, serotype) %>%
    dplyr::summarise(
      ipd_expected = sum(ipd_expected),
      .groups = "drop"
    ) %>%
    dplyr::group_by(draw) %>%
    dplyr::mutate(
      ipd_frac = ipd_expected / sum(ipd_expected)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(serotype) %>%
    dplyr::summarise(
      frac = mean(ipd_frac),
      frac_sd = stats::sd(ipd_frac),
      frac_lower = stats::quantile(ipd_frac, 0.025),
      frac_upper = stats::quantile(ipd_frac, 0.975),

      ipd_mean = mean(ipd_expected),
      ipd_expected_sd = stats::sd(ipd_expected),
      ipd_lower = stats::quantile(ipd_expected, 0.025),
      ipd_upper = stats::quantile(ipd_expected, 0.975),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      vaccine_type = dplyr::if_else(
        serotype %in% vaccine,
        "VT",
        "NVT"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(frac))

  # 8. Validate IPD compositions ----------------------------------------

  yearly_check <- ipd_yearly_summary %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      total = sum(frac),
      .groups = "drop"
    )

  stopifnot(
    all(abs(yearly_check$total - 1) < 1e-8),
    abs(sum(ipd_cumulative_summary$frac) - 1) < 1e-8
  )

  # Check that cumulative and annual expected burdens agree
  stopifnot(
    isTRUE(all.equal(
      sum(ipd_yearly_summary$ipd_mean),
      sum(ipd_cumulative_summary$ipd_mean)
    ))
  )

  # 9. Return outputs ---------------------------------------------------

  list(
    cumulative = ipd_cumulative_summary,
    yearly = ipd_yearly_summary
  )
}