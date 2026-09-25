# in this script i do things to explore whether i should be computing the Gg values 
# not from the raw metadata but accounting for the population: what fraction of the population
# with group x has the gene? if a gene is IPD-associated, or carriage-associated, Gg will 
# be a bit biased. 

# this script makes a reasonable case for this not being a problem

is_ipd <- meta_clean$collection_type == "IPD"
g <- meta_clean$group # just for short


lev  <- sort(unique(g))
ni   <- table(factor(g[is_ipd],  levels = lev))
nc   <- table(factor(g[!is_ipd], levels = lev))
both <- lev[ni >= 10 & nc >= 10]
length(both) # groups occurring in both IDP and carriage 

# cmopare the g entries for groups that are in both, but split between ipd and carriage 
gi <- rowsum(G[is_ipd, ],  g[is_ipd]);  gi <- gi / as.numeric(table(g[is_ipd])[rownames(gi)])
gc <- rowsum(G[!is_ipd, ], g[!is_ipd]); gc <- gc / as.numeric(table(g[!is_ipd])[rownames(gc)])

d <- gi[both, ] - gc[both, ]              # IPD minus carriage, per group and locus
summary(as.vector(d))

colMeans(d)[order(abs(colMeans(d)), decreasing = TRUE)][1:20]   # most divergent loci
# annoying readout though

# 
ni_b <- as.numeric(ni[both]); nc_b <- as.numeric(nc[both])
p  <- (gi[both, ] * ni_b + gc[both, ] * nc_b) / (ni_b + nc_b)   # pooled prevalence i guess that is 
# the prob of seeing the locus if you have the raw combined population 
se <- sqrt(p * (1 - p) * (1/ni_b + 1/nc_b)) # standard error you would expect in the naively-combined population
zd <- d / se # this is a two-proportion z test 
# zd should have approx standard normal distribution 

sd(zd[is.finite(zd) & se > 0])     # ~1 = entirely sampling noise; >1 = real differences

# however, note that our loci don't move independently and low values can be weird in this so: 
ok <- (p * ni_b >= 5) & ((1 - p) * ni_b >= 5) & (p * nc_b >= 5) & ((1 - p) * nc_b >= 5)
sd(zd[ok & is.finite(zd)]) # still around 1 
# suggests overall there is no bias. HOWEVER Gg might be biased for a few loci 

mz <- colMeans(zd, na.rm = TRUE)                    # per locus, averaged over groups
head(data.frame(locus = names(sort(abs(mz), decreasing = TRUE)),
                  mean_z    = round(sort(abs(mz), decreasing = TRUE), 2),
                mean_diff = round(colMeans(d)[names(sort(abs(mz), decreasing = TRUE))], 3)), 15)


# not sure those were the tests i need to do. my ai buddy suggests a non-statistical comparison 

Gg_car <- Gg
Gg_car[both, ] <- gc[both, ]                  # carriage-only values for those 43 groups

e_now <- as.numeric(x0 %*% Gg)
e_car <- as.numeric(x0 %*% Gg_car)

sum(x0[both])                                 # population weight covered by the substitution
summary(abs(e_car - e_now))                   # the bias this would correct
summary(abs(m$post2011 - m$pre2009))          # the signal you're modelling
# it does make a difference comparable to the frequencies themselves. BAH. 

# here's  a permutation-based check : 

e_now <- as.numeric(x0 %*% Gg)
# permute the IPD/carriage labels within each group. 
# This preserves group sizes and mixture proportions and destroys
# only the association, so it gives  the difference we'd see from sampling noise alone:
perm_diff <- replicate(20, {
  lab   <- ave(is_ipd, g, FUN = sample)            # shuffle within group
  gc_p  <- rowsum(G[!lab, ], g[!lab])
  gc_p  <- gc_p / as.numeric(table(g[!lab])[rownames(gc_p)])
  Gg_p  <- Gg; Gg_p[both, ] <- gc_p[both, ]
  abs(as.numeric(x0 %*% Gg_p) - e_now)
})

median(perm_diff)                       # noise-only
median(abs(e_car - e_now))              # observed
# not much. i'm going to not go into this further and I will leave G alone 

# i'll also save this as a what-about- script 
