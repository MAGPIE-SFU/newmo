# this takes meta_clean
# and reweights the counts of the groups, to adjust for the fact that the
# data have invasiveness and carriage samples from kids, adults and unknown-aged people.

# with the reweighting done, we now have the corrected locus frequencies,
# and so can make the delta statistics that we need to reorder the columns of G.
# (1) xt has the xt information (how much of each group is there at each time)

# nickstsckids = read.xlsx("../newmo_model/nick-inv-stsc-children.xlsx") # from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009389#sec027
# i did not use this one - i think the estimates were from the US?
# note that there were also 'gamma' values which were about post-vax invasiveness.
# i did not use these.
# but it might be good to do that .

#' Reweight carriage composition and define model groups
#'
#' Combines child carriage counts with IPD-derived carriage estimates (child,
#' adult and unknown-age IPD, divided by invasiveness) into a synthetic child
#' carriage population per year, defines the model groups from the
#' pre-vaccine composition, and builds x(t).
#'
#' @param meta_clean Cleaned metadata, e.g. from [clean_metadata()].
#' @param inv Invasiveness lookup, e.g. from [build_invasiveness_lookup()].
#' @param tauA Weight for adult IPD, should be very small.
#' @param tauU Weight for unknown-age IPD (proportion of the dataset that is of
#'   unknown age).
#' @param k Dirichlet prior strength: how many child carriage swabs the
#'   IPD-based estimate is worth.
#' @param prevax_years Pre-vaccine years used to define the groups.
#' @param min_sero_frac NVT serotypes below this share of pre-vaccine synthetic
#'   carriage go to "rare_NVT".
#' @param min_sg_frac Serotype-GPSC pairs below this share go to
#'   "<serotype>_other".
#'
#' @return A list with `group_lookup`, `meta_clean_groups` (meta_clean with a
#'   group column) and `final_yearly` (the synthetic
#'   carriage composition per serotype-GPSC pair, for checks).
#'
#' @importFrom dplyr %>% count filter left_join group_by mutate ungroup select
#'   transmute full_join if_else summarise bind_rows distinct across coalesce
#'   case_when rename
#' @importFrom tidyr pivot_wider replace_na complete nesting
#' @export
reweight_carriage <- function(meta_clean, inv,
                              tauA, tauU, k,
                              prevax_years, min_sero_frac, min_sg_frac) {

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

  child_v_lookup <- inv %>% transmute(serotype, v = v_child)
  adult_v_lookup <- inv %>% transmute(serotype, v = v_adult)


  child_inv <- yearly_counts %>%
      filter(age_group == "child", ipd_n > 0) %>%
      left_join(child_v_lookup, by = "serotype") %>%
      filter(!is.na(v)) %>%
      group_by(year) %>%
      mutate(        r_C = (ipd_n) / v, # carriage cases you'd expect from this number of IPD
          p_C_inv = r_C / sum(r_C)    ) %>%  # proportion group g in child carriage based on this IPD sample
      ungroup() %>%
      select( year, serotype, gpsc,
          ipd_n, carriage_n, v, r_C, p_C_inv    )

  adult_inv <- yearly_counts %>%
      filter(age_group == "adult", ipd_n > 0) %>%
      left_join(adult_v_lookup, by = "serotype") %>%
      filter(!is.na(v)) %>%
      group_by(year) %>%
      mutate( r_A = (ipd_n ) / v, # carriage cases you'd expect from this number of IPD; have removed alpha - not needed
          p_A_inv = r_A / sum(r_A)    ) %>%  # proportion group g in adult carriage based on this IPD sample
      ungroup() %>%
      select( year, serotype, gpsc,
          ipd_n, carriage_n, v, r_A, p_A_inv    )

  # unknowns: i treat them as if they are kids from the point of view of invasiveness. but we don't know.
  #  ??? might want to average 0.8*kid invasiveness + 0.2* adult invasiveness later
  unknown_proxy <- yearly_counts %>%
      filter(age_group == "unknown", ipd_n > 0) %>%
      left_join(child_v_lookup, by = "serotype") %>%
      filter(!is.na(v)) %>%
      group_by(year) %>%
      mutate(
          r_U = (ipd_n ) / v,
          p_U_proxy = r_U / sum(r_U)
      ) %>%
      ungroup() %>%
      select(year, serotype, gpsc,
          ipd_n, carriage_n, v, r_U, p_U_proxy)

  # i want to also include the carriage data we do have, such as it is..
  # here are the raw carriage fractions by group in kids

  child_carriage_counts <- yearly_counts %>%
    filter(age_group == "child") %>%
    select(year, serotype, gpsc,  carriage_n)

  # this is the data frame that is going to be turned into my final sythn population:
  components_yearly <- child_inv %>%
    select(year, serotype, gpsc,  p_C_inv) %>%
    full_join(
      adult_inv %>%
        select(year, serotype, gpsc,  p_A_inv),
      by = c("year", "serotype", "gpsc")
    ) %>%
    full_join(
      unknown_proxy %>%
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

  final_yearly <- components_yearly %>%
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
      omega_t    = N_car / (N_car + k),                           # weight on carriage (0 if none)
      p_final    = (carriage_n + k * p_C_aug) / (N_car + k),      # Dirichlet posterior mean
      p_final_sd = sqrt(p_final * (1 - p_final) / (N_car + k + 1)) # posterior SD (optional)
        ) %>%
    ungroup()

  # "Child carriage counts were combined with the IPD-derived carriage estimate
  # using a Dirichlet–multinomial model, with the IPD-derived estimate as the prior
  # mean and prior strength k (in equivalent carriage samples)."

  # now we define groups. Before, I was using a cutoff for how many isolates there were. But what if there
  # were isolates of a low-inv type? also i don't want rareNVT to be conditioned on being rare after vaccination
  # because hte NFDS model will not be able to fit that - see notes on the github from Oct 2 (CC).

  pcv7_types  <- VACCINE_TYPES$PCV7
  pcv13_extra <- setdiff(VACCINE_TYPES$PCV13, VACCINE_TYPES$PCV7)

  pre_yrs <- intersect(prevax_years, unique(final_yearly$year))

  # 1. pre-vaccine synthetic carriage share of each serotype-GPSC pair (mean over pre-vaccine years)
  pre_sg <- final_yearly %>%
    filter(year %in% pre_yrs) %>%
    group_by(year, serotype, gpsc) %>%
    summarise(p = sum(p_final), .groups = "drop") %>%
    complete(year = pre_yrs, nesting(serotype, gpsc), fill = list(p = 0)) %>%
    group_by(serotype, gpsc) %>%
    summarise(p_pre_sg = mean(p), .groups = "drop") %>%
    group_by(serotype) %>%
    mutate(p_pre_sero = sum(p_pre_sg)) %>%
    ungroup()

  # 2. a group for every serotype-GPSC pair seen at any time
  sero_tab <- distinct(pre_sg, serotype, p_pre_sero)
  sg_tab   <- select(pre_sg, serotype, gpsc, p_pre_sg)

  group_lookup <- bind_rows(distinct(meta_clean, serotype, gpsc),
                            distinct(final_yearly, serotype, gpsc)) %>%
    distinct() %>%
    left_join(sero_tab, by = "serotype") %>%
    left_join(sg_tab,   by = c("serotype", "gpsc")) %>%
    mutate(across(c(p_pre_sero, p_pre_sg), ~ coalesce(.x, 0)),      # not seen pre-vaccine -> 0
           vaccine_type = case_when(serotype %in% pcv7_types  ~ "PCV7",
                                    serotype %in% pcv13_extra ~ "PCV13",
                                    TRUE                      ~ "NVT"),
           group = case_when(
             vaccine_type == "NVT" & p_pre_sero < min_sero_frac ~ "rare_NVT",
             is.na(gpsc) | p_pre_sg < min_sg_frac             ~ paste(serotype, "other", sep = "_"),
             TRUE                                             ~ paste(serotype, gpsc, sep = "_")))

  # 3. apply the new groups, then build xt exactly as before
  final_yearly <- final_yearly %>%
    left_join(select(group_lookup, serotype, gpsc, group), by = c("serotype", "gpsc"))
  meta_clean <- meta_clean %>%
    left_join(select(group_lookup, serotype, gpsc, group), by = c("serotype", "gpsc"))

  stopifnot(!any(is.na(final_yearly$group)), !any(is.na(meta_clean$group)))

  list(
    group_lookup      = group_lookup,
    meta_clean_groups = meta_clean,
    final_yearly      = final_yearly
  )
}
