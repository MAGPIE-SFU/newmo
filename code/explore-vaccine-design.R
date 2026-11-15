# explore-vaccine-design

# going to look at how i might design a vaccine for SA given i have this nice NFDS model

# would be interesting to see whether it really depends much on the fitted parameter values

# first, how invasive are the strains and how prevalent are they?
library(dplyr); library(ggplot2); library(patchwork)


common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

xt <- read_csv(
  common$XT_PATH,
  show_col_types = FALSE
)


PCV7 = common$VACCINE_TYPES$PCV7
PCV13 = common$VACCINE_TYPES$PCV13
PCV15 = common$VACCINE_TYPES$PCV15
PCV20 = common$VACCINE_TYPES$PCV20
PCV21 = common$VACCINE_TYPES$PCV21
PCV23 = common$VACCINE_TYPES$PCV23 



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
setdiff(sero_car$serotype, sero_inv$serotype)   # rare_NVT, NTs - no invasiveness
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

# new 
compodat <- xt %>%
  filter(year %in% 2000:2009) %>%        # e.g. pre-vaccine only
  group_by(year, serotype) %>%
  summarise(frac = sum(frac), .groups = "drop") %>%   # over groups
  group_by(serotype) %>%
  summarise(frac = mean(frac), .groups = "drop")  %>%
  inner_join(sero_inv, by = "serotype") %>%
  mutate(ipdfrac = v_child*frac)%>% # does not sum to 1, but comparable from simulation to simulation
  mutate(ipdshare = ipdfrac / sum(ipdfrac)) # normalized, sums to 1

compodat <- compodat %>% mutate(vax = case_when(serotype %in% PCV7 ~ "PCV7", 
                                                serotype %in% PCV13~  "PCV13",
                                                serotype %in% PCV15 ~ "PCV15",
                                                serotype %in% PCV20~ "PCV20", # not a childhood vaccine
                                                serotype %in% PCV23 ~ "PCV23",
                                                .default = "NVT"))
vax_cols <- c(
  PCV7  = "#440154",   # dark purple
  PCV13 = "#414487",
  PCV15 = "#2a788e",
  PCV20 = "#22a884",
  PCV23 = "#7ad151",   # bright green
  NVT   = "grey70"
) 



V_FLOOR <- 10^floor(log10(min(dat$v_child)))   # 1e-05 here

# invasiveness, for reference, new version 
p_inv <- ggplot(compodat, aes(y = serotype, colour = vax)) +
  geom_segment(aes(x = V_FLOOR, xend = v_child, yend = serotype),
               linewidth = 4, lineend = "butt") +
  scale_colour_manual(values = vax_cols, name = NULL) +
  scale_x_log10(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(V_FLOOR, max(compodat$v_child))) +
  labs(x = "invasiveness (children)", y = NULL) +
  theme_bw()

p_ipdfrac <- ggplot(compodat, aes( y = serotype, colour = vax)) +
  geom_segment(aes(x = 0, xend = ipdfrac, yend = serotype),
               linewidth = 4, lineend = "butt") +
  scale_x_log10(expand = expansion(mult = c(0, 0.02))) +
    scale_colour_manual(values = vax_cols, name = NULL) +
   labs(x = "fraction of inferred pre-vaccine ipd", y = NULL) +
  theme_bw()
p_inv + p_ipdfrac + plot_layout(guides = "collect")
# end new 
# these plots illustrate the relationship between what STs are in what 
# vaccines, and their (1) child invasiveness and (2) contribution to IPD . 
# potentially useful supplemental figure for illustrating why we choose the candidates we do
# in fact it would be even better to make another colour for being a new candidate ST


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


p_inv + p_carpre+ p_change + plot_layout(guides = "collect")

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

candidates = setdiff(summary$serotype[1:27] , common$VACCINE_TYPES$PCV7) # could compute this directly from inv, might be more reliable than making my summary
current   <- common$VACCINE_TYPES$PCV7
remaining <- candidates
trace     <- tibble(step = integer(), added = character(),
                    ipd = double(), ipd_drop = double())
targetsize <- 20 # valency
baseIPD <- score_vaccine(modelpars, current, inv, group_lookup)

while (length(current) < targetsize && length(remaining) > 0) {

  scores <- vapply(remaining,
                   function(s) score_vaccine(modelpars, c(current, s), inv, group_lookup),
                   FUN.VALUE = 1)

  best <- remaining[which.min(scores)]

  trace <- add_row(trace,
                   step     = length(current) - 7 + 1, # 7 because it's the length of pcv7..
                   added    = best,
                   ipd      = min(scores),
                   ipd_drop = baseIPD - min(scores))

  current   <- c(current, best)
  remaining <- setdiff(remaining, best)
  base_ipd  <- min(scores)

  message(sprintf("step %d: added %s, IPD %.1f",
                  nrow(trace), best, base_ipd))
}
trace
#  the IPD does not drop much as you add the next few
# note that we do better than PCV13; can cmopare to our ipdfrac15 vaccine for example



# now let's do an exhaustive search
library(tictoc)
tic(); score_vaccine(modelpars, current, inv, group_lookup); toc() # 0.03s per evaluation, can do 32 per second
# so an exhaustive search of the top 25 serotypes, but forced to include PCV7 (so 18 to choose from),
# and we search for the remaining 8 to add for a 15-valent vaccine, is choose(18,8) which will take about 1350s
# or about 30 mins. seems reasonable.


library(future.apply)

plan(multisession, workers = parallel::detectCores(logical = FALSE) - 2)   # 16
candidates = setdiff(summary$serotype[1:27] , common$VACCINE_TYPES$PCV7) # could compute this directly from inv, might be more reliable than making my summary
combos <- combn(candidates, 8, simplify = FALSE)   # combos[[1]] is a character vector of things added to PCV7
# the vaccine for each row of the resulting stuff is going to be c(PCV7, combos[[id]])

# compute all the scores
scores <- future_vapply(combos,
                        function(s) score_vaccine(modelpars, c(PCV7, s), inv, group_lookup),
                        FUN.VALUE=1)

res8 <- tibble(
  id    = seq_along(combos),
  added = vapply(combos, paste, character(1), collapse = ","),
  ipd   = scores
) %>%
  arrange(ipd) %>%
  mutate(rank = row_number(),
         ipd_rel = ipd / ipd[1])

write_csv(res8, "exhaustive_15valent.csv")
plan(sequential)       # release the workers
