library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate) # may not be needed 
library(openxlsx) # read.xlxs


# this script takes meta_clean 
# and reweights the counts of the groups, to adjust for the fact that the 
# data have invasiveness and carriage samples from kids, adults and unknown-aged people. 

# with the reweighting done, we now have the corrected locus frequencies, 
# and so can make the delta statistics that we need to reorder the columns of G. 
# this script creates information for 2 files (and saves them): 
# (1) xt_grouped.csv has the xt information (how much of each group is there at each time)
# (2) locus_freqs_deltas.csv has the locus frequency information for each year and the "delta stastic"


# nickstsckids = read.xlsx("../newmo_model/nick-inv-stsc-children.xlsx") # from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009389#sec027
# i did not use this one - i think the estimates were from the US? 
# note that there were also 'gamma' values which were about post-vax invasiveness.
# i did not use these.
# but it might be good to do that . 

# pseudocounts and weights 
tauA <- 0.01 # weight for adult IPD, should be v small 
tauU <- (272/(3865+1386)) # prop of the dataset that are unknown age - small; REDO if new dataset obvs

# sept 22: omega was the weight for carriage compared to IPD. I set it to 0.27 because that's the fraction of the data that's from carriage 
# note that we have 439 adult carriage samples that aren't contributing 
# here even though adult IPD is seen as reflective of the carriage in the population
# should we use kids only throughout , or adults only throughout? or this weird mix? 
# seems like carriage is directly what we want and IPD is a biased sample - up-weight carriage? 
# setting tauA very small so adult IPD contributes little. LL's paper said other studies only measured kids. 
# Sept 24: did some digging about better ways to do this. 
# will use a way with a better statistical interpretation (works out very similarly) 
k <- 780   # Dirichlet prior strength- interp is how many child carriage swabs the IPD-based
# estimate is worth. Carriage weight in a year = N / (N + k);
# k = 780 gives ~0.27 at N ~ 290 (matches omega is about 0.27 - one par instead of 2 and can depend on the number of carriage samples).


######################################################
# Begin: count all the stuff and figure out what portion of the carriage population each group 
# contributes in each year. this will eventually be x(t) for the NFDS model. 
######################################################

meta_clean <- read_csv("data/SouthAfrica/meta_clean.csv")
inv <- read_csv("data/Invasiveness/invasiveness_lookup.csv") # see make_invasiveness_lookup.R script

yearly_counts_long <- meta_clean %>%
    count( year, age_group, collection_type,
        serotype, gpsc, group,  name = "n" )

hist(filter(yearly_counts_long, year==2007)$n,breaks = 50) # in a given year most are 1 (!)
# which probably means i need to figure out which types to remove
# and also means i should have alpha << 1 (alpha is a parameter to avoid problems with 0s later. so is beta)

# make it wider by splitting IPD and carriage to cols
yearly_counts <- yearly_counts_long %>%
    pivot_wider(
        names_from = collection_type,
        values_from = n,
        values_fill = 0  ) %>%  rename( ipd_n = IPD,carriage_n = carriage )

# connect it to the supplementary tables from Nick's paper via the inv data frame 
missing <- setdiff(unique(meta_clean$serotype), inv$serotype)
if (length(missing)) stop("No invasiveness for: ", paste(missing, collapse = ", "),
                          " -- re-run make_invasiveness_lookup.R")

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
    select( year, serotype, gpsc, group,
        ipd_n, carriage_n, v, r_C, p_C_inv    )

adult_inv <- yearly_counts %>%
    filter(age_group == "adult", ipd_n > 0) %>%
    left_join(adult_v_lookup, by = "serotype") %>%
    filter(!is.na(v)) %>%
    group_by(year) %>%
    mutate( r_A = (ipd_n ) / v, # carriage cases you'd expect from this number of IPD; have removed alpha - not needed  
        p_A_inv = r_A / sum(r_A)    ) %>%  # proportion group g in adult carriage based on this IPD sample
    ungroup() %>%
    select( year, serotype, gpsc, group,
        ipd_n, carriage_n, v, r_A, p_A_inv    )

# unknowns: i treat them as if they are kids from the point of view of invasiveness. but we don't know. See notes / sanity check at the bottom of this script. 
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
    select(year, serotype, gpsc, group,
        ipd_n, carriage_n, v, r_U, p_U_proxy)

# i want to also include the carriage data we do have, such as it is.. 
# here are the raw carriage fractions by group in kids 
# beta helps avoid 0s


# 
child_carriage_counts <- yearly_counts %>%
  filter(age_group == "child") %>%
  select(year, serotype, gpsc, group, carriage_n)

# this is the data frame that is going to be turned into my final sythn population: 
components_yearly <- child_inv %>%
  select(year, serotype, gpsc, group, p_C_inv) %>%
  full_join(
    adult_inv %>%
      select(year, serotype, gpsc, group, p_A_inv),
    by = c("year", "serotype", "gpsc", "group")
  ) %>%
  full_join(
    unknown_proxy %>%
      select(year, serotype, gpsc, group, p_U_proxy),
    by = c("year", "serotype", "gpsc", "group")
  ) %>%
  full_join(
    child_carriage_counts,
    by = c("year", "serotype", "gpsc", "group")
  ) %>%
  mutate(
    p_C_inv = replace_na(p_C_inv, 0),
    p_A_inv = replace_na(p_A_inv, 0),
    p_U_proxy = replace_na(p_U_proxy, 0),
    carriage_n = replace_na(carriage_n, 0)
  ) %>%
   ungroup()

# sanity check: do the things add to 1 where they are supposed to? 
child_inv %>%
    group_by(year) %>%
    summarise(sum_p_C_inv = sum(p_C_inv))

adult_inv %>%
    group_by(year) %>%
    summarise(sum_p_A_inv = sum(p_A_inv))

unknown_proxy %>%
    group_by(year) %>%
    summarise(sum_p_U_proxy = sum(p_U_proxy))


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
# I do not know the best values of the pars tauA, tauU and omega

# Leonie's paper notes that the NFDS model datasets to date only had data from children (!) 
# suggests that tauA (the adult contribution) should be basically 0 to be in keeping with those studies 
# since we have it, i will use it a tiny bit, i suppose. 
# but we might want to simplify our model by removing adults entirely before model-building

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
# mean and prior strength k (in equivalent carriage samples)." k is defined above. 

head(final_yearly)

# sanity check -- the sums should be 1 
final_yearly %>%
    group_by(year) %>%
    summarise(sum_p_C_aug = sum(p_C_aug),
        sum_p_final = sum(p_final)  ) # yes they all sum to 1. 


ggplot(final_yearly, aes(x=p_final, fill=year))+geom_histogram() +facet_wrap(~as.factor(year))+xlim(c(0, 0.1))




###########################################################################
# make x(t) and save . 
###########################################################################

xt <- final_yearly %>%
  group_by(year, group) %>%
  summarise(frac = sum(p_final), a0 = first(N_car) + k, .groups = "drop") %>%
  complete(year, group, fill = list(frac = 0)) %>%
  group_by(year) %>% mutate(a0 = max(a0, na.rm = TRUE)) %>% ungroup() %>%
  mutate(frac_sd = sqrt(frac * (1 - frac) / (a0 + 1)),
         count   = frac * 1e5) %>%
  left_join(distinct(meta_clean, group, vaccine_type), by = "group") %>%
  select(year, group, vaccine_type, frac, frac_sd, count) %>%
  arrange(year, group)

# some checks 
xt %>% group_by(year) %>% summarise(s = sum(frac)) %>% pull(s) %>% range()   # should be 1 1

# (that's the smallest and largest total fraction in a year - should both be 1, all years are ok) 
sum(is.na(xt$vaccine_type)) # nothing has no vaccine type listed, good,  should be 0

xt %>% group_by(year, vaccine_type) %>% summarise(f = sum(frac), .groups = "drop") %>%
  pivot_wider(names_from = vaccine_type, values_from = f) # VT should fall after 2009 and 2011, mostly does 

glimpse(xt)

write_csv(xt, "data/SouthAfrica/xt_grouped.csv")


################################################################################
# Aside: compare the unknown age people to kids and to adults; let's see which they look more like. i bet kids .
################################################################################
group_var <- "serotype"

comp_overall_sero <- meta_clean %>%
  count(collection_type, age_group, .data[[group_var]], name = "n") %>%
  group_by(collection_type, age_group) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  rename(grouping = all_of(group_var))

comp_wide_sero <- comp_overall_sero %>%
  select(collection_type, age_group, grouping, p) %>%
  pivot_wider(
    names_from = age_group,
    values_from = p,
    values_fill = 0
  )

comp_wide_sero %>%
  group_by(collection_type) %>%
  summarise(
    TV_unknown_child = 0.5 * sum(abs(unknown - child)),
    TV_unknown_adult = 0.5 * sum(abs(unknown - adult)),
    cor_unknown_child = cor(unknown, child),
    cor_unknown_adult = cor(unknown, adult)
  )
ggplot(comp_wide_sero, aes(x = child, y = unknown)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~collection_type) +
  theme_bw() +
  labs(
    x = "Child composition",
    y = "Unknown-age composition"
  )
ggplot(comp_wide_sero, aes(x = adult, y = unknown)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~collection_type) +
  theme_bw() +
  labs(
    x = "Adult composition",
    y = "Unknown-age composition"
  )
# the unknowns look more like kids than adults in their ST composition. 
