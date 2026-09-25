# locus frequency exploration and plotting 

# this should not really be in the Ggrouping script, so i have moved it here 
# but if, in rebuilding such a model, the sanity checks here are not sensible,
# that indicates a problem that should be managed in one of: 
# - the choice of accessory genes to include (and merge, etc) (Tanu) 
# - the choice of gruops to keep or not 
# - the choice of k in the Dirichlet thing in reweighting 


###########################################################################
# locus frequencies across the years, all plotted against 2000 
###########################################################################


# ---- PLOT - ---
# reorg to long form so i can make a plot with frequency in the year 2000 on the x
# axis and frequency in other years on the y axis, for 2001 to 2014 
plotdat <- as.data.frame(locus_by_year[, 1:15] ) %>%
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


# CONCLUSION: overall good stability in the locus frequencies. 
# in this weighting, locus frequencies are very stable until about 2013, 2014. 
# they probably change more during vaccination than before
# NOTE that there were post-vax different invasiveness values in Nick's paper 
# (which i think are the gamma values in Nick's paper) ?
# but i don't want to revisit that now.

# vaccination: PCV7 in 2009, PCV13 in 2011 (See Bellman paper figure 1). 



# - aside 
# do we get the delta info from just pre-vax frequency comparison? 
tmp <- locus_by_year %>% 
  transmute(
    pre2005 = rowMeans(across(num_range("t", 2000:2004))),
    f2005.2008 = rowMeans(across(num_range("t", 2005:2008))),
    delta = (post2011 - pre2009)^2/(1-pre2009*(1-pre2009)),
    delPRE = (f2005.2008 - pre2005)^2/(1-pre2005*(1-pre2005))
  )
ggplot(tmp, aes(x = sqrt(delta), y = sqrt(delPRE)  )) + geom_point(alpha = 0.1) 
# interesting. No, we don't get 
# the delta information from pre-vaccine frequency comparisons. but there is something



# and i want to know more about z and delta, now that i have all this other information in the locus by year data 

