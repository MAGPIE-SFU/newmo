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

# eventually maybe we should do an inverse variance weighting or something. 
# but for now i'll use Nick's inv. estimates from his paper in a series of adjustments for age and carriage vs ipd to get synthetic carriage population

nickkids = read.xlsx("data/Invasiveness/nick-inv-sts-children.xlsx")
nickadults = read.xlsx("data/Invasiveness/nick-inv-sts-adults.xlsx")

# nickstsckids = read.xlsx("../newmo_model/nick-inv-stsc-children.xlsx") # from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009389#sec027
# i did not use this one - i think the estimates were from the US? 
# note that there were also 'gamma' values which were about post-vax invasiveness.
# i did not use these.
# but it might be good to do that . 

# pseudocounts 
alpha <- 0.05          # pseudocount in r_g = (d_g + alpha)/v_g -- need to check that we aren't adding a ton of wrongness with alpha too high
omega <- 0.10          # carriage weight in the final blend
beta <- 0.05           # another pseudo count, in child carriage 


######################################################
# Begin: count all the stuff and figure out what portion of the carriage population each group 
# contributes in each year. this will eventually be x(t) for the NFDS model. 
######################################################

meta_clean = read_csv("data/SouthAfrica/meta_clean.csv")

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

# connect it to the supplementary tables from Nick's paper 
child_v_lookup <- nickkids %>%
    transmute(  serotype = Serotype,
        v = invasiveness  )

adult_v_lookup <- nickadults %>%
    transmute( serotype = Serotype,
        v = invasiveness  )

child_inv <- yearly_counts %>%
    filter(age_group == "child") %>%
    left_join(child_v_lookup, by = "serotype") %>%
    filter(!is.na(v)) %>%
    group_by(year) %>%
    mutate(        r_C = (ipd_n + alpha) / v, # carriage cases you'd expect from this number of IPD
        p_C_inv = r_C / sum(r_C)    ) %>%  # proportion group g in child carriage based on this IPD sample
    ungroup() %>%
    select( year, serotype, gpsc, group,
        ipd_n, carriage_n, v, r_C, p_C_inv    )

adult_inv <- yearly_counts %>%
    filter(age_group == "adult") %>%
    left_join(adult_v_lookup, by = "serotype") %>%
    filter(!is.na(v)) %>%
    group_by(year) %>%
    mutate( r_A = (ipd_n + alpha) / v, # carriage cases you'd expect from this number of IPD 
        p_A_inv = r_A / sum(r_A)    ) %>%  # proportion group g in adult carriage based on this IPD sample
    ungroup() %>%
    select( year, serotype, gpsc, group,
        ipd_n, carriage_n, v, r_A, p_A_inv    )

# unknowns: i treat them as if they are kids from the point of view of invasiveness. but we don't know. See notes / sanity check at the bottom of this script. 
#  ??? might want to average 0.8*kid invasiveness + 0.2* adult invasiveness later  
unknown_proxy <- yearly_counts %>%
    filter(age_group == "unknown") %>%
    left_join(child_v_lookup, by = "serotype") %>%
    filter(!is.na(v)) %>%
    group_by(year) %>%
    mutate(
        r_U = (ipd_n + alpha) / v,
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
  group_by(year) %>%
  mutate(
    total_carriage = sum(carriage_n),
    q_car_child = if_else(
      total_carriage > 0,
      (carriage_n + beta) / sum(carriage_n + beta), # fraction of this group in the carriage data 
      NA_real_
    )
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
tauA <- 0.01
tauU <- (272/(3865+1386)) # prop of the dataset that are unknown age - small 

omega <- 0.27 # omega is the contribution of our carriage dataset 
# i suppose perhaps it should be more than 0.1 (where we hvae carriage..) 
# sept 22: setting to 0.27 because that's the fraction of the data that's from carriage 
# note that we have 439 adult carriage samples that aren't contributing 
# here even though adult IPD is seen as reflective of the carriage in the population
# should we use kids only throughout , or adults only throughout? or this weird mix? 
# seems like carriage is directly what we want and IPD is a biased sample - up-weight carriage? 


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
    
    has_child_carriage = sum(carriage_n) > 0,
    p_final = if_else(
      has_child_carriage,
      (1 - omega) * p_C_aug + omega * q_car_child,
      p_C_aug
    )
  ) %>%
  ungroup()


head(final_yearly)

# sanity check -- the sums should be 1 
final_yearly %>%
    group_by(year) %>%
    summarise(sum_p_C_aug = sum(p_C_aug),
        sum_p_final = sum(p_final)  ) # yes they all sum to 1. 


ggplot(final_yearly, aes(x=p_final, fill=year))+geom_histogram() +facet_wrap(~as.factor(year))+xlim(c(0, 0.1))




###########################################################################
# make x(t) and start exploring locus and ST frequencies over time 
###########################################################################

# NOTE - sept 22 changed this file, will push to git, but i did NOT
# re-save the reweighted output or recreate G from Tanu's file so from 
# here, this has to be completed 
# and it would be good to have an order in which to run all this to recreate with different values 
# of the fixed parameters 

final_yearly = final_yearly %>% mutate( syncount = p_final*100000)
# x(t) is now just final_yearly$syncount for year t 

xt = final_yearly %>% select(year, serotype, gpsc,group, p_final, syncount) %>%
    rename(frac = p_final, count = syncount) 

# we want the order to match our grouped G matrix, so we can set up our dynamics 
# without having to worry about reordering 

Gg = read.csv("data/SouthAfrica/Ggrouped.csv", row.names = 1, check.names = FALSE)

xt <- xt %>% mutate(group = factor(group, levels = rownames(Gg))) %>% 
    complete(year, group, fill = list(frac = 0, count = 0)) %>% 
    arrange(year, group)

# sanity check? 
all(xt$group[1:50] == rownames(Gg)[1:50]) # TRUE 

write_csv(xt, "data/SouthAfrica/xt_grouped.csv")


###########################################################################
# locus frequencies across the years, all plotted against 2000 
###########################################################################


years <- sort(unique(final_yearly$year))

locus_by_year <- sapply(years, function(y) {
    p_df <- final_yearly %>%
        filter(year == y) %>%
        select(group, p_final)
    p_vec <- setNames(p_df$p_final, p_df$group)[rownames(G_mean)]
    p_vec[is.na(p_vec)] <- 0
    as.numeric(t(p_vec) %*% as.matrix(Gg))
})

rownames(locus_by_year) <- colnames(Gg)
colnames(locus_by_year) <- paste("t",years, sep="")
locus_by_year = as.data.frame(locus_by_year)
plot(locus_by_year$t2000, locus_by_year$t2014)

# ---- PLOT - ---
# reorg to long form so i can make a plot with frequency in the year 2000 on the x
# axis and frequency in other years on the y axis, for 2001 to 2014 
plotdat <- as.data.frame(locus_by_year) %>%
    mutate(locus = rownames(locus_by_year)) %>%
    pivot_longer(
        cols = -locus,
        names_to = "year",
        values_to = "freq"
    ) %>%
    left_join(
        data.frame(
            locus = rownames(locus_by_year),
            freq_2000 = locus_by_year[, "t2000"]
        ),
        by = "locus"
    ) %>%
    filter(year != "t2000")

ggplot(plotdat, aes(x = freq_2000, y = freq)) +
    geom_point(alpha = 0.5, size = 0.8,color="blue") +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~ year) +
    theme_bw() +
    labs(
        x = "Locus frequency in 2000",
        y = "Locus frequency in year T"
    )


# CONCLUSION: in this weighting, locus frequencies are very stable until about 2013, 2014. 
# i wonder if this is because i did not account for post-PCV changes in the  invasiveness 
# (which i think are the gamma values in Nick's paper) ? but i don't want to revisit that now. i want to use G
# and get on with making a model and comparing it to data. 

# vaccination: PCV7 in 2009, PCV13 in 2011 (See Bellman paper figure 1). 

# Corander et al defined Delta as (f - e)^2 / (1 - e(1-e)) 
# where e is the frequency before vaccination and f is the frequency after vaccination 
# (all pre- and during vaccination samples for e, all post-vaccination samples for f). 

# here, to do approximately the same, i would take the pre-2009 and post-2011 
# mean frequencies of these 

locus_by_year <- locus_by_year %>%
    mutate(
       pre2009  = rowMeans(across(num_range("t", 2000:2008))),
       post2011 = rowMeans(across(num_range("t", 2012:2014))),
    )

plot(locus_by_year$pre2009, locus_by_year$post2011) # sensible 


# - aside 
# do we get the delta info from just pre-vax frequency comparison? 
tmp <- locus_by_year %>% 
    transmute(
        pre2005 = rowMeans(across(num_range("t", 2000:2004))),
        f2005.2008 = rowMeans(across(num_range("t", 2005:2008))),
        delta = (post2011 - pre2009)^2/(1-pre2009*(1-pre2009)),
        delPRE = (f2005.2008 - pre2005)^2/(1-pre2005*(1-pre2005))
    )

plot(tmp$delta, tmp$delPRE  ) 
# interesting. No, we don't get 
# the delta information from pre-vaccine frequency comparisons. 

locus_by_year <- locus_by_year %>% 
    mutate( 
        delta = (post2011 - pre2009)^2/(1-pre2009*(1-pre2009))
            )

hist(locus_by_year$delta,breaks = 50)
plot(locus_by_year$pre2009, locus_by_year$delta) 
# no clear relationship 
# between frequency pre-vaccination and delta. 

# ok, let's proceed to use these deltas. 

# want to order Gg's columns from lowest to highest delta 
# so i might as well order the loci this way here too 
locsinorder <- rownames(locus_by_year)[order(locus_by_year$delta)]

locus_by_year <- locus_by_year[ order(locus_by_year$delta) , ] 

write_csv(locus_by_year, "data/SouthAfrica/locus_freqs_deltas.csv")





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