# File:    sensitivity_test.R

# Assess sensitivity of the NFDS model to its input parameters by perturbing
#  inputs near the optimized parameter settings


# Preliminaries
# ==============================================================================

library(ggplot2)
library(purrr)

# load the simulate() function for solving the NFDS model
source("simulate_ode.R")

# load in some common variable settings
source("common.R")

# load the fitting functions
source("parameter_fitting_functions.R")

# load in fits
opt_pars <- readRDS("../data/SouthAfrica/all_best_fitted_parameters.rds")

# default theme for the plots
clean_theme <- theme(panel.background = element_rect(fill="white", colour="black"),
      panel.grid.major = element_line("grey85"),
      panel.grid.minor = element_line('grey95'))


# some relevant parameters for running the test and generating figures
N_GROUPS <- 25     # number of groups to plot
QTILE <- 0.95      # summary quantile for the ribbons
N_SIMS <- 1000     # number of simulated parameter settings to draw (for each loss function)
START_YEAR <- 2008 # year prior to start of model evaluation (i.e., fitted trajectories start at START_YEAR+1)


# Optimal trajectories 
# ==============================================================================

# get the model
sa_model <- get_samodel()

top_groups <- data.frame(IC = sa_model$x0) |>
  rownames_to_column(var="group") |>
  mutate(frac = IC / sum(IC)) |>
  arrange(desc(frac)) |>
  slice_head(n=N_GROUPS)

# process results for each of the three optimal parameter settings
out1 <- run_model(as.list(opt_pars[1,]), T) |>
  rownames_to_column(var="group") |>
  # filter(group %in% top_groups$group) |>
  pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
  bind_cols(loss = opt_pars[1,"loss"])

out2 <- run_model(as.list(opt_pars[2,]), T) |>
  rownames_to_column(var="group") |>
  # filter(group %in% top_groups$group) |>
  pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
  bind_cols(loss = opt_pars[2,"loss"])

out3 <- run_model(as.list(opt_pars[3,]), T) |>
  rownames_to_column(var="group") |>
  # filter(ST %in% top_STs$ST) |>
  pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
  bind_cols(loss = opt_pars[3,"loss"])

# plot trajectories for the three optimal fits
bind_rows(out1, out2, out3) |>
  filter(group %in% top_groups$group) |>
  ggplot(aes(x=year, y=value, colour=loss)) +
  geom_point() +
  geom_line() +
  facet_wrap(~group, ncol=5, scales="free_y") +
  clean_theme +
  theme(legend.position="bottom")



# Perturbations
# ==============================================================================

# For each parameter simulate n_sims draws from a normal distribution with mean
#  equal to the optimized values and standard deviation loosely based on the 
#  loss function profiles

sd_pars <- data.frame(par = c("log10v1", "log10sigma", "log10migr", "portion", "log10whigh", "log10wfrac"),
                      SD  = c( 0.25,      0.5,          1,           0.1,       0.25,         2),
                      # bounds on the sampling scale: portion in [0,1], and
                      #  log10wfrac <= 0 so that wlow <= whigh
                      lo  = c(-Inf,      -Inf,         -Inf,         0,        -Inf,         -Inf),
                      hi  = c( Inf,       Inf,          Inf,         1,         Inf,          0))

# floor for wfrac so that fits with wlow = 0 have a finite mean on the log10
#  scale (log10(0) = -Inf breaks the truncated normal draws)
wfrac_min <- 1e-8

par_lookup <- lapply(1:3, function(i) {as.data.frame(opt_pars[i,]) |>
    mutate(log10v1 = log10(v1),
           log10sigma = log10(sigma),
           log10migr = log10(migr),
           log10whigh = log10(whigh),
           log10wfrac = log10(pmax(wlow/whigh, wfrac_min))) |>
    select(log10v1, log10sigma, log10migr, log10whigh, log10wfrac, portion) |>
    pivot_longer(everything(), names_to="par") |>
    left_join(sd_pars, by="par")
  }
)

# draw from a normal truncated to [lo, hi] via the inverse CDF
sim_pars <- Vectorize(function(lookup, par_name){
  row <- lookup[lookup$par == par_name, ]
  p <- runif(1, pnorm(row$lo, row$value, row$SD), pnorm(row$hi, row$value, row$SD))
  qnorm(p, row$value, row$SD)},
  vectorize.args="par_name")

par_df <- lapply(1:3, function(j) {
  lapply(1:N_SIMS, function(i){
    data.frame(par = sd_pars$par) |>
      mutate(val = sim_pars(par_lookup[[j]], par)) |>
      bind_cols(rep=i)
  }) |> 
    list_rbind() |>
    bind_cols(loss=names(loss_functions)[j])
}) |>
  list_rbind()

sim_par_sets <- par_df |>
  pivot_wider(names_from="par", values_from = "val") |>
  mutate(sigma = 10^log10sigma,
         v1 = 10^log10v1,
         migr = 10^log10migr,
         whigh = 10^log10whigh,
         wfrac = 10^log10wfrac,
         .keep = "unused")

stopifnot(all(between(sim_par_sets$portion, 0, 1)),
          all(sim_par_sets$wfrac <= 1))

# run the model for each perturbed parameter set, keeping the rep/loss labels
#  alongside the trajectories so the results can be summarised and plotted

# PSOCK cluster (works on Windows); leave one core free for the OS
n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)

# each worker is a fresh R session, so it needs the same working directory and
#  the same sourced model code as this one
parallel::clusterCall(cl, setwd, getwd())
parallel::clusterEvalQ(cl, {
  source("simulate_ode.R")
  source("common.R")
  source("parameter_fitting_functions.R")
  NULL
})
parallel::clusterExport(cl, c("sim_par_sets", "top_groups"))

run_one <- function(k) {
  pars_k <- sim_par_sets[k, ]
  # a bad parameter draw shouldn't kill the whole batch; drop it and record why
  tryCatch(
    run_model(as.list(select(pars_k, -rep, -loss)), T) |>
      rownames_to_column(var="group") |>
      # filter(group %in% top_groups$ST) |>
      pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
      bind_cols(select(pars_k, rep, loss)),
    error = function(e) {
      data.frame(rep = pars_k$rep, loss = pars_k$loss, error = conditionMessage(e))
    }
  )
}

# load-balanced in chunks: ODE solve times vary a lot across parameter draws,
#  so hand out small batches rather than one fixed block per core
sim_res <- parallel::parLapplyLB(cl, seq_len(nrow(sim_par_sets)), run_one,
                                 chunk.size = 25)
parallel::stopCluster(cl)

failed <- keep(sim_res, \(x) "error" %in% names(x)) |> list_rbind()
if (NROW(failed) > 0) warning(nrow(failed), " parameter sets failed; see `failed`")

sim_out <- discard(sim_res, \(x) "error" %in% names(x)) |>
  list_rbind()

# plot perturbed trajectories for each loss function
sim_out |>
  filter(group %in% top_groups$group) |>
  group_by(group, year, loss) |>
  summarise(L = quantile(value, probs=(1-QTILE)/2),
            med = quantile(value, probs=0.5),
            U = quantile(value, probs=1-(1-QTILE)/2),
            .groups = "drop") |>
  mutate(year = START_YEAR + year) |>
  ggplot(aes(x=year, colour=loss)) +
  geom_ribbon(aes(ymin=L, ymax=U, fill=loss), alpha=0.25, colour=NA) +
  geom_line(aes(y=med), linewidth=1) +
  facet_wrap(~group, ncol=5, scales="free_y") +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function") +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  clean_theme +
  theme(legend.position="bottom")



# Combine the above plot with xt and the optimal fits
# ==============================================================================

# data
xt <- read_csv(common$XT_PATH, col_types=cols())

# fitted trajectories
opt_fit <- bind_rows(out1, out2, out3) |>
  filter(group %in% top_groups$group) |>
  mutate(year = START_YEAR + year)

# simulated trajectory summaries
sim_summ <- sim_out |>
  filter(group %in% top_groups$group) |>
  group_by(group, year, loss) |>
  summarise(L = quantile(value, probs=(1-QTILE)/2),
            med = quantile(value, probs=0.5),
            U = quantile(value, probs=1-(1-QTILE)/2),
            .groups = "drop") |>
  mutate(year = START_YEAR + year)

xt |>
  filter(group %in% top_groups$group) |>
  ggplot(aes(x=year)) +
  geom_point(aes(y=frac)) +
  geom_segment(aes(y=frac-frac_sd*1.96, yend=frac+frac_sd*1.96)) +
  geom_line(data=opt_fit, aes(y=value, colour=loss), linewidth=1) +
  geom_ribbon(data=sim_summ, aes(ymin=L, ymax=U, fill=loss), alpha=0.2) +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function",
       title=paste0("Top ",N_GROUPS, " ST-GPSC groups. Ribbons show ", QTILE*100, "% range.")) +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  facet_wrap(~group, scales="free_y") +
  clean_theme +
  theme(legend.position="bottom")
ggsave(paste0(here::here("plots/"), "sensitivity_top",N_GROUPS,"_ribbon",QTILE,".png"), width=9, height=7, units="in", dpi=700)
