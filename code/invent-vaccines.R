# invent vaccines based on a greedy approach, and on an exhaustive search

# prelims
library(dplyr)
library(patchwork)
library(future.apply)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

source("code/simulate_ode.R")
source("code/explore-vaccine-fns.R")
source("code/parameter_generator.R")


# Load information
samodel <- readRDS(common$SA_MODEL_PATH)
xt <- read_csv(
  common$XT_PATH,
  show_col_types = FALSE
)
nyears <- common$NYEARS

inv <- read.csv("data/SouthAfrica/processed-inputs/invasiveness_lookup_interval.csv")
group_lookup <- read.csv("data/SouthAfrica/processed-inputs/grouplookup.csv")

PCV7 = common$VACCINE_TYPES$PCV7
PCV13 = common$VACCINE_TYPES$PCV13
PCV15 = common$VACCINE_TYPES$PCV15
PCV21 = common$VACCINE_TYPES$PCV21
PCV23 = common$VACCINE_TYPES$PCV23

# rds_params_fit_path <- common$ALL_BEST_FITS_PATH
rds_params_fit_path <- common$OVERALL_BEST_FITS_PATH

all_fits <- readRDS(rds_params_fit_path)

# ===========================================================================
# set up modelpars: relevant for methods 2 and 3, not 1
# ===========================================================================

# best overall fit
best_fit_model_pars <- readRDS(common$OVERALL_BEST_FITS_PATH)

modelpars <- parameter_generator(
  how_numerical = best_fit_model_pars$loss,
  optimizer = best_fit_model_pars$method,
  fit_id = best_fit_model_pars$fit_id,
  path = common$OVERALL_BEST_FITS_PATH,
  vaccine_type1 = "PCV7", vaccine_type2 = "PCV13",
  how_gene_order = "default",
  use_switch = TRUE
)

modelpars$tswitch <- 0 # we won't use a switch
modelpars$gene_order <- 1:length(samodel$el) # logitdiff order since that is how we made Gnfds

# ---- method 1: use just the estimated contribution to the IPD population in our inferred dataset
sero_inv <- inv %>%
  select(serotype, v_child) %>%
  filter(!is.na(v_child))

#data on the composition of the inferred carriage population before vaccination, including 2009
compodat <- xt %>%
  filter(year %in% 2000:2009) %>%        # e.g. pre-vaccine only
  group_by(year, serotype) %>%
  summarise(frac = sum(frac), .groups = "drop") %>%   # over groups
  group_by(serotype) %>%
  summarise(frac = mean(frac), .groups = "drop")  %>%
  inner_join(sero_inv, by = "serotype") %>%
  mutate(ipdfrac = v_child*frac)%>% # does not sum to 1, but comparable from simulation to simulation
  mutate(ipdshare = ipdfrac / sum(ipdfrac)) # normalized, sums to 1

# the top 25 serotypes, as a character vector, in order of the fraction of IPD they are causing in our
# inferred pre-vaccine population

# what i would want here is the 13-valent vaccine that has the top contributors to IPD
top13 <- compodat %>%
  slice_max(ipdfrac, n = 13) %>%
  pull(serotype) %>%
  as.character() # as it happens this does include PCV7

vax_pre_ipdfrac = top13 # the vaccine you would choose in method 1.
# we could also add sequentially to pcv7 adding one of these at a time but let's keep 13 valent vaccines for now


# what about the top 13 invasive strains? (not fraction contribs)?
top13car <- compodat %>%
  slice_max(v_child, n = 13) %>%
  pull(serotype) %>%
  as.character()

direct = tibble(vaccine = compress_vstring(vax_pre_ipdfrac),
                ipd = score_vaccine(modelpars, vax_pre_ipdfrac, inv, group_lookup, years=6:10),
                method = "ipdfrac")
direct = rbind(direct,
               tibble(vaccine = compress_vstring(top13car),
                      ipd = score_vaccine(modelpars, top13car, inv, group_lookup, years=6:10),
                      method = "mostinvasive"))

#####################
# aside (for now) 
#####################
# what are the top 25 most invasive strains and are pcv15, 21 adn 23 subsets of these? 
# if not, what rank do you have to get to to get all the ones in these pcvs? 
top25 <- compodat %>%
  slice_max(v_child, n = 25) %>%
  pull(serotype) %>%
  as.character()
setdiff(top25, PCV23) # in top25 but not pcv23
unique (c(setdiff(top25, PCV23), setdiff(PCV23, top25))) # 

# will define a new candidate set 
compodat <- compodat %>% mutate(vax = case_when(serotype %in% PCV7 ~ "PCV7", 
                                  serotype %in% PCV13~  "PCV13",
                                  serotype %in% PCV15 ~ "PCV15",
                                  serotype %in% PCV21~ "PCV21",
                                  serotype %in% PCV23 ~ "PCV23",
                                  .default = "NVT"))

ggplot(compodat, aes(y = serotype, colour = vax)) +
  geom_segment(aes(x = V_FLOOR, xend = v_child, yend = serotype),
               linewidth = 4, lineend = "butt") +
 # scale_colour_manual(values = vax_cols, name = NULL) +
  scale_x_log10(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(V_FLOOR, max(dat$v_child))) +
  labs(x = "invasiveness (children)", y = NULL) +
  theme_bw()

#####################
# direct
#####################
# ---- method 2: greedy .
# explanation: of the top 25 invasive serotypes (NOT according to fraction in carriage, just
# according to their estimated invasiveness in children), choose the best one
# to add to PCV7. Then choose the best one to add to c(PCV, firstone). Etc until
# we get to a 13-valent vaccine. Could rank by frac IPD but didn't . Probably would
# get very similar results

candidates <- compodat %>%
  filter(!serotype %in% PCV7) %>%
  slice_max(v_child, n = 18) %>%
  pull(serotype) %>%
  as.character()

current   <- common$VACCINE_TYPES$PCV7
remaining <- candidates
trace     <- tibble(step = integer(), vaccine = character(),
                    ipd = double())
targetsize <- 13 # valency
baseIPD <- score_vaccine(modelpars, current, inv, group_lookup)

tmp=PCV7
while (length(current) < targetsize && length(remaining) > 0) {

  scores <- vapply(remaining,
                   function(s) score_vaccine(modelpars, c(current, s), inv, group_lookup, years = 6:10),
                   FUN.VALUE = 1)

  best <- remaining[which.min(scores)]

  trace <- add_row(trace,
                   step     = length(current) - 7 + 1, # 7 because it's the length of pcv7..
                   vaccine    = compress_vstring(c(tmp,best)),
                   ipd      = min(scores))
  tmp = c(tmp, best)
  current   <- c(current, best)
  remaining <- setdiff(remaining, best)
  base_ipd  <- min(scores)

  message(sprintf("step %d: added %s, IPD %.1f",
                  nrow(trace), best, base_ipd))
}
# trace
# note that we do better than PCV13 by the end (a little)

# ---- method 3: exhaustive search

plan(multisession, workers = parallel::detectCores(logical = FALSE) - 2)   # 16

candidates <- compodat %>%
  filter(!serotype %in% PCV7) %>%
  slice_max(v_child, n = 18) %>%
  pull(serotype) %>%
  as.character() # same starting candidates as in method 2

# all possible 18 choose 6 combinations of the candidates
combos <- combn(candidates, 6, simplify = FALSE)   # combos[[1]] is a character vector of things added to PCV7

# for each one, add those to PCV7 and compute the score
scores <- future_vapply(combos,
                        function(s) score_vaccine(modelpars, c(PCV7, s), inv, group_lookup, years = 6:10),
                        FUN.VALUE=1)

# collect the results
exhaustive_results <- tibble(
  id    = seq_along(combos),
  vaccine = vapply(combos,
                   function(s) paste(c(PCV7, s), collapse = ","),
                   FUN.VALUE = character(1)),
  ipd   = scores
) %>%
  arrange(ipd)

plan(sequential)       # release the workers

# join all the results together
trace$method= "greedy"
exhaustive_results$method = "exhaustive"

all_vaccines <- bind_rows(
  trace               %>% select(ipd, vaccine, method),
  direct              %>% select(ipd, vaccine, method),
  exhaustive_results  %>% select(ipd, vaccine, method)
)

# and save (no long with commas but i'll stick to rds anyway)
saveRDS(all_vaccines, file = "data/SouthAfrica/processed-inputs/vaccines.rds")
