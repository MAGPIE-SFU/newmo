# Generation of xt_grouped, incorporating uncertainty about Nick's progression rates

# this takes meta_clean
# and reweights the counts of the groups, to adjust for the fact that the
# data have invasiveness and carriage samples from kids, adults and unknown-aged people.
# The groups themselves are defined in reweight_carriage(); here invasiveness is
# drawn n_reps times from its log-normal distribution and x(t) is averaged over the draws.

#' Grouped x(t) with invasiveness uncertainty
#'
#' Repeats the carriage reweighting `n_reps` times, each time with one
#' invasiveness draw per serotype (from the log-normal in the invasiveness
#' lookup), builds x(t) for the given groups, and combines the repetitions.
#'
#' @param meta_clean Cleaned metadata, e.g. from [clean_metadata()].
#' @param inv Invasiveness lookup, e.g. from [build_invasiveness_lookup()].
#' @param group_lookup Group definitions, e.g. `reweight_carriage()$group_lookup`.
#' @param tauA Weight for adult IPD, should be very small.
#' @param tauU Weight for unknown-age IPD.
#' @param k Dirichlet prior strength, in equivalent child carriage swabs.
#' @param lambda Fraction of the prior mean that is uniform.
#' @param n_reps Number of repetitions for the invasiveness uncertainty.
#' @param seed Random seed, set at the start as in the original script.
#'
#' @return x(t) with columns year, group, serotype, vaccine_type, frac, frac_sd.
#'
#' @importFrom dplyr %>% count rename transmute filter mutate select inner_join
#'   group_by ungroup full_join if_else left_join summarise n_distinct first n
#'   bind_rows any_of
#' @importFrom tidyr pivot_wider replace_na complete
#' @importFrom stats rlnorm var
#' @export
reweight_carriage_uncertainty <- function(meta_clean, inv, group_lookup,
                                          tauA, tauU, k,
                                          lambda = 0.01, n_reps = 1000, seed = 1) {

  # set the seed
  set.seed(seed)

  # Processing the invasive data
  # ==============================================================================

  ######################################################
  # Begin: count all the stuff and figure out what portion of the carriage population each group
  # contributes in each year. this will eventually be x(t) for the NFDS model.
  ######################################################

  yearly_counts_long <- meta_clean %>%
    count( year, age_group, collection_type,
           serotype, gpsc,   name = "n" )

  # make it wider by splitting IPD and carriage to cols
  yearly_counts <- yearly_counts_long %>%
    pivot_wider(
      names_from = collection_type,
      values_from = n,
      values_fill = 0  ) %>%  rename( ipd_n = IPD,carriage_n = carriage )

  # connect it to the supplementary tables from Nick's paper via the inv data frame
  missing <- setdiff(unique(meta_clean$serotype), inv$serotype)
  if (length(missing)) stop("No invasiveness for: ", paste(missing, collapse = ", "),
                            " -- rebuild the invasiveness lookup")

  # NEW   - take the log normal pars from the table before rather than moment matching
  # to Nick's ranges, because in the invasiveness lookup
  # we did "shrinkage" for the log means
  child_v_lookup <- inv %>%
      transmute(serotype, logmean = v_child_logmean, logsd = v_child_logsd) %>%
      filter(!is.na(logmean), !is.na(logsd))

  adult_v_lookup <- inv %>%
      transmute(serotype, logmean = v_adult_logmean, logsd = v_adult_logsd) %>%
      filter(!is.na(logmean), !is.na(logsd))

  # one draw per serotype, reused across all years within a replicate
  draw_v <- function(lk) lk %>% mutate(v = rlnorm(n(), logmean, logsd)) %>% select(serotype, v)


  get_child_inv <- function(vdraw) {
    yearly_counts %>%
      filter(age_group == "child", ipd_n > 0) %>%
      inner_join(vdraw, by = "serotype") %>%
      group_by(year) %>% # CC: previously each st got a draw of v for each year? instead we want ONE draw for the ST, not one for the st for each SC and year (!)
      mutate( r_C = (ipd_n) / v, # carriage cases you'd expect from this number of IPD
             p_C_inv = r_C / sum(r_C)) %>%  # proportion group g in child carriage based on this IPD sample
      ungroup() %>%
      select(year, serotype, gpsc,
             ipd_n, carriage_n, v, r_C, p_C_inv)
  }

  get_adult_inv <- function(vdraw) {
    yearly_counts %>%
      filter(age_group == "adult", ipd_n > 0) %>%
      inner_join(vdraw, by = "serotype") %>%
      group_by(year) %>%
      mutate(r_A = (ipd_n ) / v, # carriage cases you'd expect from this number of IPD; have removed alpha - not needed
             p_A_inv = r_A / sum(r_A)) %>%  # proportion group g in adult carriage based on this IPD sample
      ungroup() %>%
      select(year, serotype, gpsc,
             ipd_n, carriage_n, v, r_A, p_A_inv)
  }


  # unknowns: i treat them as if they are kids from the point of view of invasiveness. but we don't know.
  #  ??? might want to average 0.8*kid invasiveness + 0.2* adult invasiveness later
  get_unknown_proxy <- function(vdraw) {
    yearly_counts %>%
      filter(age_group == "unknown", ipd_n > 0) %>%
      inner_join(vdraw, by = "serotype") %>%
      group_by(year) %>%
      mutate(r_U = (ipd_n ) / v,
             p_U_proxy = r_U / sum(r_U)) %>%
      ungroup() %>%
      select(year, serotype, gpsc,
             ipd_n, carriage_n, v, r_U, p_U_proxy)
  }


  # i want to also include the carriage data we do have, such as it is..
  # here are the raw carriage fractions by group in kids

  child_carriage_counts <- yearly_counts %>%
    filter(age_group == "child") %>%
    select(year, serotype, gpsc,  carriage_n)

  # this is the data frame that is going to be turned into my final sythn population:
  components_yearly <- lapply(1:n_reps, function(x) {
      vc <- draw_v(child_v_lookup)      # child AND unknown share this draw
      va <- draw_v(adult_v_lookup)

      get_child_inv(vc) %>%
    select(year, serotype, gpsc,  p_C_inv) %>%
    full_join(
      get_adult_inv(va) %>%
        select(year, serotype, gpsc,  p_A_inv),
      by = c("year", "serotype", "gpsc")
    ) %>%
    full_join(
      get_unknown_proxy(vc) %>%
        select(year, serotype, gpsc,  p_U_proxy),
      by = c("year", "serotype", "gpsc")
    ) %>%
    full_join(
      child_carriage_counts,
      by = c("year", "serotype", "gpsc")
    ) %>%
    mutate(
      p_C_inv = replace_na(p_C_inv, 0),
      p_A_inv = replace_na(p_A_inv, 0),
      p_U_proxy = replace_na(p_U_proxy, 0),
      carriage_n = replace_na(carriage_n, 0)
    ) %>%
    ungroup()
  })

  ######################################################################
  # create the final weighted composition accounting for all the various things
  ######################################################################

  # now, recognizing that not all the years have any carriage data, i want to create my final weighted synthetic carriage
  # population for all of my years

  # first i need to choose tauU and tauA, and create
  # p_C_aug = (1-tauU-tauA)*p_C_inv + tauA*p_A_inv + tau_U*p_U_proxy

  # p_C is the portion of this group in the (synthetic) child carriage population estimated from inv disease data
  # it is a balance of what you'd estimate from the invasive disease data in children (using Nick's invasiveness),
  # the inv disease in adults, and the inv disease from the unknown age people.

  # then i can combine that with my carriage sample for the years when i have carriage

  final_yearly <- lapply(components_yearly, function(x) {
    x %>%
      group_by(year) %>%
      mutate(
        has_A = sum(p_A_inv) > 0,
        has_U = sum(p_U_proxy) > 0,

        tauA_eff = if_else(has_A, tauA, 0),
        tauU_eff = if_else(has_U, tauU, 0),

        tau_total = (1 - tauA - tauU) + tauA_eff + tauU_eff,

        wC = (1 - tauA - tauU) / tau_total,
        wA = tauA_eff / tau_total,
        wU = tauU_eff / tau_total,

        p_C_aug = wC * p_C_inv + wA * p_A_inv + wU * p_U_proxy,
        N_car      = sum(carriage_n),                               # child carriage swabs this year
        omega_t    = N_car / (N_car + k)                           # weight on carriage (0 if none)
    # CC: removed the p_final and p_final_sd lines here, pseudocount was pulling too far towards uniform by x 84 in the denominator
       ) %>%
      ungroup()
  })

  # "Child carriage counts were combined with the IPD-derived carriage estimate
  # using a Dirichlet–multinomial model, with the IPD-derived estimate as the prior
  # mean and prior strength k (in equivalent carriage samples)."

  # 3. apply the groups, then build xt
  final_yearly <- lapply(final_yearly, function(x){
    x %>%
      left_join(select(group_lookup, serotype, gpsc, group), by = c("serotype", "gpsc"))
  })

  meta_clean <- meta_clean %>%
      select(-any_of("group")) %>% # in case group was already in there i guess?
      left_join(select(group_lookup, serotype, gpsc, group), by = c("serotype", "gpsc"))

  stopifnot(!any(sapply(final_yearly, function(x) any(is.na(x$group)))), !any(is.na(meta_clean$group)))

  ###########################################################################
  # make x(t). but keep st in xt
  ###########################################################################
  group_sero <- group_lookup %>%
    group_by(group) %>%
    summarise(serotype     = if (n_distinct(serotype) == 1) first(serotype) else first(group),
              vaccine_type = first(vaccine_type),
              .groups = "drop")

  all_groups <- sort(unique(group_lookup$group))
  all_years  <- sort(unique(final_yearly[[1]]$year))

  # uses lambda to smooth the prior ps, not the counts
  xt <- lapply(final_yearly, function(x) {
      x %>%
          group_by(year, group) %>%
          summarise(p_tilde_g  = sum(p_C_aug),       # sums to 1 over groups
                    carriage_g = sum(carriage_n),    # sums to N_car
                    .groups = "drop") %>%
          complete(year = all_years, group = all_groups,
                   fill = list(p_tilde_g = 0, carriage_g = 0)) %>%
          group_by(year) %>%
          mutate(N_car     = sum(carriage_g),
                 p_tilde   = (1 - lambda) * p_tilde_g + lambda / n(),   # n() == 84
                 alpha     = carriage_g + k * p_tilde,
                 frac_mean = alpha / sum(alpha),                        # sum(alpha) == N_car + k
                 frac_var  = frac_mean * (1 - frac_mean) / (N_car + k + 1)) %>%
          ungroup() %>%
          left_join(group_sero, by = "group") %>%
          select(year, group, serotype, vaccine_type, frac_mean, frac_var)
  })

  # aggregate results and compute the new standard deviation
  bind_rows(xt) %>%
    group_by(year, group, serotype, vaccine_type) %>%
    summarise(frac = mean(frac_mean),
              frac_sd = sqrt(mean(frac_var) + var(frac_mean)), .groups="drop")
}


################################################################################
### plot of the estimated carriage composition with the uncertainty
################################################################################

#' Plot the most common groups in x(t) with uncertainty
#'
#' @param xt_agg x(t) from [reweight_carriage_uncertainty()].
#' @param n_top Number of most common groups to show.
#' @return A ggplot.
#' @importFrom dplyr %>% group_by summarise slice_max pull filter left_join mutate
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_line geom_point geom_hline
#'   geom_vline facet_wrap scale_x_continuous labs theme_grey theme element_text
#'   element_rect element_blank
#' @export
plot_xt_top_groups <- function(xt_agg, n_top = 20) {

  # Select the most common groups
  top_groups <- xt_agg %>%
    group_by(group) %>%
    summarise(mean_frac = mean(frac), .groups = "drop") %>%
    slice_max(mean_frac, n = n_top, with_ties = FALSE) %>%
    pull(group)

  # Calculate pre-vaccine mean for each group
  prevax_means <- xt_agg %>%
    filter(year < 2010) %>%
    group_by(group) %>%
    summarise(
      prevax_mean = mean(frac),
      .groups = "drop"
    )

  # Prepare plotting data
  plot_data <- xt_agg %>%
    filter(group %in% top_groups) %>%
    left_join(prevax_means, by = "group") %>%
    mutate(
      lower = pmax(0, frac - 1.96 * frac_sd),
      upper = pmin(1, frac + 1.96 * frac_sd),
      group = factor(group, levels = top_groups)
    )

  # Plot
  ggplot(plot_data, aes(x = year, y = frac)) +

    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0,
      colour = "grey65",
      linewidth = 0.5
    ) +

    geom_line(
      colour = "grey70",
      linewidth = 0.4
    ) +

    geom_point(
      colour = "black",
      size = 1.6
    ) +

    geom_hline(
      aes(yintercept = prevax_mean),
      linetype = "dashed",
      colour = "grey35",
      linewidth = 0.5
    ) +

    geom_vline(
      xintercept = 2010,
      linetype = "dotted",
      linewidth = 0.5
    ) +

    facet_wrap(~ group, scales = "free_y", ncol = 5) +

    scale_x_continuous(
      breaks = seq(2000, 2020, 5)
    ) +

    labs(
      title = paste("The", n_top, "most common ST-SC groups"),
      subtitle = "Estimated carriage composition with 95% uncertainty intervals",
      x = NULL,
      y = "Frequency",
      caption = "Points: estimated carriage frequency. Grey bars: approximate 95% intervals. Dashed: pre-vaccine mean."
    ) +

    theme_grey(base_size = 11) +

    theme(
      strip.text = element_text(
        colour = "royalblue",
        face = "bold",
        size = 10
      ),
      strip.background = element_rect(fill = "grey85"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
}
