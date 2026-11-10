# locus frequency exploration and plotting 

# this should not really be in the Ggrouping script, so i have moved it here 
# but if, in rebuilding such a model, the sanity checks here are not sensible,
# that indicates a problem that should be managed in one of: 
# - the choice of accessory genes to include (and merge, etc) (Tanu) 
# - the choice of gruops to keep or not 
# - the choice of k in the Dirichlet thing in reweighting 


# read in some of the necessary stuff 
xt <- read_csv("data/SouthAfrica/processed-inputs/xt_grouped.csv", show_col_types = FALSE) 
G <- as.matrix(read.csv("data/SouthAfrica/processed-inputs/Gnfds.csv", row.names = 1))
locus_by_year <- read_csv("data/SouthAfrica/processed-inputs/locus_freqs_deltas.csv")

rownames(locus_by_year) = locus_by_year$locus

###########################################################################
# locus frequencies across the years, all plotted against 2000 
###########################################################################


# ---- PLOT - ---
# reorg to long form so i can make a plot with frequency in the year 2000 on the x
# axis and frequency in other years on the y axis, for 2001 to 2014 
plotdat <- as.data.frame(locus_by_year[, 1:15] ) %>%
  mutate(locus = rownames(locus_by_year)) %>%
  pivot_longer(
    cols = 1:15,
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
glimpse(plotdat)


ggplot(plotdat, aes(x = t2000, y = freq)) +
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
# the delta information from pre-vaccine frequency comparisons.
# but there is something
# i wonder if we get any of the others like z or logitdiff 





# and i want to know more about z and delta, now that i have all this other information in the locus by year data 

ggplot(locus_by_year, aes(x=z, y=delta))+geom_point(alpha=0.3)
ggplot(locus_by_year, aes(x=abs(z), y=delta, colour = pre2009))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()


# alternatives log(f/e)  and logit(f) - logit(e) (both | | of those obvs) 
# for less overall tendency to choose high vs low frequency genes 
# logit <- function(x) log(x/(1-x)) 
# locus_by_year <- locus_by_year %>% 
#  mutate( lograt = abs(log(pre2009/post2011)) ) 

ggplot(filter(locus_by_year, in_model), aes(x=abs(pre2009), y=lograt, color = logitdiff))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()
# log ratio is good too 

ggplot(filter(locus_by_year, in_model == TRUE) , aes(x=abs(pre2009), y=logitdiff))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()
# this has very low frequency dependence. Note that the high-freq genes will have to be pretty stable so
# it requires a huge lineage turnover to shift them. (though they don't have the highest delta)
# choosing the lowest lograt genes for the NFDS model will, even with the 0.95 cutoff,
# lead us to quite a few high-frequency genes. whereas logitdiff is not going to bias towards low or high -freq genes

ggplot(filter(locus_by_year, in_model == TRUE), 
       aes(x=logitdiff, y = lograt, color = pre2009)) + 
    geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()

# wow, there is a population of genes for which the log ratio and logitdiff do not agree
# and these are the higher-frequency genes 

ggplot(filter(locus_by_year, in_model == TRUE), 
       aes(x=delta, y = z, color = pre2009)) + 
    geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()
# delta and z are very similar 

# note also: consider whether 0 or 1 are so different. should we want our measure to be 
# symmetric under the transformation p -> 1-p, such that we think NFDS also applies to 
# the lack of a gene? The log ratio does *not* have this property but delta, z and the logit do. 

ggplot(filter(locus_by_year, in_model == TRUE) , aes(x=delta, y=logitr, color = pre2009))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()

# how do the rankings actually compare? 
m <- locus_by_year[locus_by_year$in_model, ]
n <- 200                                      # or whatever the model uses
length(intersect(head(m$locus[order(m$delta)],  n),
                 head(m$locus[order(m$logitdiff)], n)))
cor(m$delta, m$logitdiff, method = "spearman")

# where they are different what are the frequencies? 
swapped <- setdiff(head(m$locus[order(m$logitdiff)], 200),
                   head(m$locus[order(m$delta)],  200))
summary(m$pre2009[m$locus %in% swapped])
summary(m$pre2009)
# logitdiff will choose more intermediate-frequency loci than delta will. which i think is what we want. 
# it is symmetric under the swap p to 1-p which treats absence of a locus as a trait too 
# and its slope against pre-vaccine frequency is v small - frequency-independent. 
# therefore it seems a better choice than delta , but not vastly 




  # some slopes: shows that delta goes up with frequency, z and logitr are approx stable and lograt goes down. 
 m <- locus_by_year[locus_by_year$in_model, ]
 sapply(c("delta", "z", "logitdiff", "lograt"), function(v)
   coef(lm(log(abs(m[[v]])) ~ log(m$pre2009)))[2])
# logitdiff now the only one with a near-zero slope 
 
# --- now i want to explore whether i can get a good NFDS result out of pre-vaccine stability measures 
 
 # so let's compute these measures over the 9 years of pre-vaccine data 
 # i can compute the variance over time for each gene and use that 
 # and that should be the variance in the logitdiff 
 # but i can compute the variance in the other things too : delta, z, whatever . 
 
 # already have some stuff in here like varianc in the frequency 
 ggplot(filter(locus_by_year, in_model), aes(x=obs_pre, y=logitdiff))+geom_point(alpha=0.5) + 
      scale_x_log10()+ scale_y_log10()
# varlogitpre is no good - interesting but could see if any of the others 
# none of the others are correlated with logitdiff either  
