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

ggplot(locus_by_year, aes(x=z, y=delta))+geom_point(alpha=0.3)
ggplot(locus_by_year, aes(x=abs(z), y=delta, colour = pre2009))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()


# alternatives log(f/e)  and logit(f) - logit(e) (both | | of those obvs) 
# for less overall tendency to choose high vs low frequency genes 
logit <- function(x) log(x/(1-x)) 
locus_by_year <- locus_by_year %>% 
  mutate( lograt = abs(log(pre2009/post2011)), 
          logitr = abs(logit(pre2009)-logit(post2011) ) ) 

ggplot(locus_by_year, aes(x=abs(pre2009), y=lograt))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()


ggplot(filter(locus_by_year, in_model == TRUE) , aes(x=abs(pre2009), y=lograt))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()
# this is still frequency dependent because of course the high-freq genes will have to be pretty stable
# it requires a huge lineage turnover to shift them. (though they don't have the highest delta)
# choosing the lowest lograt genes for the NFDS model will, even with the 0.95 cutoff,
# lead us to quite a few high-frequency genes. 


# note also: consider whether 0 or 1 are so different. should we want our measure to be 
# symmetric under the transformation p -> 1-p, such that we think NFDS also applies to 
# the lack of a gene? The log ratio does *not* have this property but delta, z and the logit do. 

# let's take a look at the logit 
ggplot(filter(locus_by_year, in_model == TRUE) , aes(x=abs(pre2009), y=logitr))+ 
geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()
# abs(logit  difference) is looking like a winner here: choosing the most stable will give us
# genes across a range of frequencies, not just low-frequency genes which i thikn will happen with 
# delta and z. no aactually delta is ok here i geuss. 
 
ggplot(filter(locus_by_year, in_model == TRUE) , aes(x=delta, y=logitr, color = pre2009))+ 
  geom_point(alpha=0.5) + scale_x_log10()+ scale_y_log10()

# how do the rankings actually compare? 
m <- locus_by_year[locus_by_year$in_model, ]
n <- 200                                      # or whatever the model uses
length(intersect(head(m$locus[order(m$delta)],  n),
                 head(m$locus[order(m$logitr)], n)))
cor(m$delta, m$logitr, method = "spearman")

# where they are different what are the frequencies? 
swapped <- setdiff(head(m$locus[order(m$logitr)], 200),
                   head(m$locus[order(m$delta)],  200))
summary(m$pre2009[m$locus %in% swapped])
summary(m$pre2009)
# logitr will choose more intermediate-frequency loci than delta will. which i think is what we want. 
# it is symmetric under the swap p to 1-p which treats absence of a locus as a trait too 
# and its slope against pre-vaccine frequency is v small - frequency-independent. 
# therefore it seems a better choice than delta , but not vastly 




  # some slopes: shows that delta goes up with frequency, z and logitr are approx stable and lograt goes down. 
 m <- locus_by_year[locus_by_year$in_model, ]
 sapply(c("delta", "z", "logitr", "lograt"), function(v)
   coef(lm(log(abs(m[[v]])) ~ log(m$pre2009)))[2])




