# explore-vaccine-design 

# going to look at how i might design a vaccine for SA given i have this nice NFDS model 

# would be interesting to see whether it really depends much on the fitted parameter values

# first, how invasive are the strains and how prevalent are they? 
library(dplyr); library(ggplot2); library(patchwork)

# --- serotype-level carriage fraction -----------------------------------
# averaged over years; change the filter for a single year or a period
sero_car <- xt_agg %>%
   filter(year %in% 2000:2008) %>%        # e.g. pre-vaccine only
  group_by(year, serotype) %>%
  summarise(frac = sum(frac), .groups = "drop") %>%   # over groups
  group_by(serotype) %>%
  summarise(frac = mean(frac), .groups = "drop")      # over years

sero_post <- xt_agg %>% 
  filter(year %in% 2013:2014) %>%        # e.g. pre-vaccine only
  group_by(year, serotype) %>%
  summarise(frac = sum(frac), .groups = "drop") %>%   # over groups
  group_by(serotype) %>%
  summarise(frac = mean(frac), .groups = "drop")      # over years

sero_inv <- inv %>%
  select(serotype, v_child) %>%
  filter(!is.na(v_child))

# mis-matches? 
setdiff(sero_car$serotype, sero_inv$serotype)   # rare_NVT, NTs, *_other labels
setdiff(sero_inv$serotype, sero_car$serotype)   # not in the population

dat <- inner_join(sero_car, sero_inv, sero_post, by = "serotype") %>%
  mutate(serotype = factor(serotype, levels = serotype[order(v_child)]))

VAX = common$VACCINE_TYPES$PCV13
dat <- inner_join(sero_car, sero_inv, by = "serotype") %>%
  mutate(vax = if_else(serotype %in% VAX, "PCV13", "NVT"),
         serotype = factor(serotype, levels = serotype[order(v_child)]))

dat <- sero_inv %>%
  left_join(sero_car,  by = "serotype") %>% rename(frac_pre  = frac) %>%
  left_join(sero_post, by = "serotype") %>% rename(frac_post = frac) %>%
  filter(!is.na(frac_pre) | !is.na(frac_post)) %>%
  mutate(across(c(frac_pre, frac_post), ~ tidyr::replace_na(.x, 0)),
         vax = if_else(serotype %in% VAX, "PCV13", "NVT"),
         serotype = factor(serotype, levels = serotype[order(v_child)]),
         ipdfrac = frac_pre*v_child)

vax_cols <- c(PCV13 = "#c0392b", NVT = "grey45")

V_FLOOR <- 10^floor(log10(min(dat$v_child)))   # 1e-05 here

# invasiveness, for reference
p_inv <- ggplot(dat, aes(y = serotype, colour = vax)) +
  geom_segment(aes(x = V_FLOOR, xend = v_child, yend = serotype),
               linewidth = 4, lineend = "butt") +
  scale_colour_manual(values = vax_cols, name = NULL) +
  scale_x_log10(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(V_FLOOR, max(dat$v_child))) +
  labs(x = "invasiveness (children)", y = NULL) +
  theme_bw()

# portion in IPD, mean, from carriage samples, note could gete this from the metadata too and compare
p_ipdfrac <- ggplot(dat, aes(x=ipdfrac,y=serotype, fill=vax))+
  geom_col() +
  scale_fill_manual(values = vax_cols, guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "fraction of inferred pre-vaccine carriage", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# portion carriage, mean over the pre-vaccine years 
p_carpre <- ggplot(dat, aes(frac_pre, serotype, fill = vax)) +
  geom_col() +
  scale_fill_manual(values = vax_cols, guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "fraction of inferred pre-vaccine carriage", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# portion carriage, post-vaccine 
p_carpost <- ggplot(dat, aes(frac_post, serotype, fill = vax)) +
  geom_col() +
  scale_fill_manual(values = vax_cols, guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "fraction of inferred post-vaccine carriage", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# really nice plot showing arrows for change in the % in carriage, log scale. 
p_change <- ggplot(dat, aes(y = serotype)) +
  geom_segment(aes(x = frac_pre, xend = frac_post, yend = serotype,
                   colour = vax),
               linewidth = 1,
               arrow = arrow(length = unit(0.07, "in"), type = "closed")) +
  geom_point(aes(x = frac_pre), colour = "grey40", size = 1.6) +
  scale_colour_manual(values = vax_cols, name = NULL) +
  scale_x_log10(labels = scales::percent) +
  labs(x = "carriage fraction, pre → post", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


p_inv + p_ipdfrac+ p_change + plot_layout(guides = "collect")

# note: what we see is that PCV13 does pretty well at choosing the top contributors to IPD. 
# there are several alternate strains we could put in if we want the top 13 contributors to IPD. 
# we can see this here: 
dat2 <- dat %>%  mutate( serotype = factor(serotype, levels = serotype[order(ipdfrac)]))
ggplot(dat2, aes(x=ipdfrac,y=serotype, fill=vax))+
  geom_col() +
  scale_fill_manual(values = vax_cols, guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "fraction of inferred pre-vaccine carriage", y = NULL) +
  theme_bw() 
# this visually shows the difference between PCV13 and the top predicted IPD contributors 
# using a mean of the pre-vaccine data (but we've gone through the inference process and out again;
# could compare to raw IPD samples instead but i'm not right now) 

ind = order(dat2$ipdfrac, decreasing = TRUE)
summary = data.frame(  serotype = as.character(dat2$serotype)[ind], VT = dat2$vax[ind], ipdfrac = dat2$ipdfrac[ind])
summary

################################################################################################
# ---- IPD-fraction-based alternatives to PCV13 for SA. or to PCVk where there are k serotypes in the vaccine
# the k-valent vaccine consisting of the top k contributors to IPD in our mean pre-vaccine population 
# is the vaccine consisting of the top k serotypes in this summary list. 

ipdfrac13 = summary$serotype[1:13]
ipdfrac15 = summary$serotype[1:15]
ipdfrac20 = summary$serotype[1:20] 


################################################################################################
# --- NFDS principles help us design a vaccine
################################################################################################ 

# i think i do this: start with PCV7. for each ST in the top 25 invasiveness list, add that to PCV7 if it isn't in there already. 
# pick the best one for minimizing IPD, call it S1 
# now start with PCV + S1, and for the ST in the remaining list, add each one to PCV7 + S1. Pick the best and call it S2. 
# continue until you have a k-valent vaccine 

# for this i will need: set up modelpars; set the vaccine; simulate; getIPD for the simulation, etc 







