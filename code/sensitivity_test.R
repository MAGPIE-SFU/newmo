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
opt_pars <- readRDS(paste0(here::here("data/SouthAfrica"), "/all_fitted_parameters_switch.rds"))

# making modifications for the initial file Lars sent me
opt_pars <- opt_pars |>
  filter(start_id == "round_06") |>
  bind_cols(whigh=1)

# default theme for the plots
clean_theme <- theme(panel.background = element_rect(fill="white", colour="black"),
      panel.grid.major = element_line("grey85"),
      panel.grid.minor = element_line('grey95'))

# pull out initial conditions
IC <- data.frame(value = get_samodel()$x0) |>
  rownames_to_column(var="group") |> 
  bind_cols(year=0)


# some relevant parameters for running the test and generating figures
N_GROUPS <- 16       # number of groups to plot
QTILE <- 0.95        # summary quantile for the ribbons
N_SIMS <- 100        # number of simulated parameter settings to draw (for each loss function)
START_YEAR <- 2007.5 # year prior to start of model evaluation (i.e., fitted trajectories start at START_YEAR+1)
BEST_LOSS <- "SSE"   # the loss function used to obtain the optimal parameters (options are SSE, weighted_SSE, and weighted_composition_distance)


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
  bind_rows(IC) |>
  bind_cols(loss = opt_pars[1,"loss"])

out2 <- run_model(as.list(opt_pars[2,]), T) |>
  rownames_to_column(var="group") |>
  # filter(group %in% top_groups$group) |>
  pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
  bind_rows(IC) |>
  bind_cols(loss = opt_pars[2,"loss"])

out3 <- run_model(as.list(opt_pars[3,]), T) |>
  rownames_to_column(var="group") |>
  # filter(ST %in% top_STs$ST) |>
  pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
  bind_rows(IC) |>
  bind_cols(loss = opt_pars[3,"loss"])

# plot trajectories for the three optimal fits
bind_rows(out1, out2, out3) |>
  filter(group %in% top_groups$group) |>
  ggplot(aes(x=year, y=value, colour=loss)) +
  geom_point() +
  geom_line() +
  facet_wrap(~group, ncol=4, scales="free_y") +
  clean_theme +
  theme(legend.position="bottom")



# Perturbations
# ==============================================================================

# For each parameter simulate n_sims draws from a normal distribution with mean
#  equal to the optimized values and standard deviation loosely based on the 
#  loss function profiles

sd_pars <- data.frame(par = c("logv1", "logsigma", "portion"),
                      SD  = c( 0.1,      0.25,          0.01),
                      # bounds on the sampling scale: portion in [0,1], and
                      lo  = c(-Inf,      -Inf,         0.01),
                      hi  = c( Inf,       Inf,         0.99))

par_lookup <- lapply(1:3, function(i) {as.data.frame(opt_pars[i,]) |>
    mutate(logv1 = log(v1),
           logsigma = log(sigma),
           portion = portion,
           .keep = "none") |>
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

# pull out parameters used in the optimal fits
opt_settings <- lapply(seq_along(par_lookup), function(i) {
  par_lookup[[i]] |>
    select(par, value) |>
    bind_cols(loss = names(loss_functions)[i])
}) |> list_rbind()

# check the range of the simulated parameters (on the simulation scale)
par_df |>
  ggplot(aes(x=val)) +
  geom_histogram(aes(fill=loss), colour="black") +
  facet_grid(loss~par, scales="free") +
  geom_vline(data=opt_settings, aes(xintercept=value), linewidth=1, linetype="dashed") +
  clean_theme

sim_par_sets <- par_df |>
  pivot_wider(names_from="par", values_from = "val") |>
  mutate(sigma = exp(logsigma),
         v1 = exp(logv1),
         .keep = "unused")

stopifnot(all(between(sim_par_sets$portion, 0, 1)))

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
parallel::clusterExport(cl, c("sim_par_sets", "top_groups", "IC"))

run_one <- function(k) {
  pars_k <- sim_par_sets[k, ]
  # a bad parameter draw shouldn't kill the whole batch; drop it and record why
  tryCatch(
    run_model(as.list(select(pars_k, -rep, -loss)), T) |>
      rownames_to_column(var="group") |>
      # filter(group %in% top_groups$ST) |>
      pivot_longer(-group, names_to="year", names_prefix="year.", names_transform = as.numeric) |>
      bind_rows(IC) |>
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

# fitted trajectories
opt_fit <- bind_rows(out1, out2, out3) |>
  filter(group %in% top_groups$group) |>
  mutate(year = START_YEAR + year)

# Plot: ribbons w/ optimal trajectories (dashed line)
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
  geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
  facet_wrap(~group, ncol=4, scales="free_y") +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function") +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  clean_theme +
  theme(legend.position="bottom")

# plot all simulated trajectories
sim_out |>
  filter(group %in% top_groups$group) |>
  mutate(year = START_YEAR + year) |>
  ggplot(aes(x=year, colour=loss)) +
  geom_line(aes(y=value, colour=loss, group=interaction(rep,loss)), alpha=0.1) +
  geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
  facet_wrap(~group, ncol=4, scales="free_y") +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function") +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  clean_theme +
  theme(legend.position="bottom")

# plot all simulated trajectories: SSE only
sim_out |>
  filter(group %in% top_groups$group,
         loss =="SSE") |>
  mutate(year = START_YEAR + year) |>
  ggplot(aes(x=year, colour=loss)) +
  geom_line(aes(y=value, colour=loss, group=interaction(rep,loss)), alpha=0.1) +
  geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
  facet_wrap(~group, ncol=4, scales="free_y") +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function") +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  clean_theme +
  theme(legend.position="bottom")
ggsave(paste0(here::here("plots"), "/sim_trajectories_SSE.pdf"), width=9, height=7, units="in", dpi=700)

# plot all simulated trajectories: weighted SSE only
sim_out |>
  filter(group %in% top_groups$group,
         loss =="weighted_SSE") |>
  mutate(year = START_YEAR + year) |>
  ggplot(aes(x=year, colour=loss)) +
  geom_line(aes(y=value, colour=loss, group=interaction(rep,loss)), alpha=0.1) +
  geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
  facet_wrap(~group, ncol=4, scales="free_y") +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function") +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  clean_theme +
  theme(legend.position="bottom")
ggsave(paste0(here::here("plots"), "/sim_trajectories_weighted-SSE.pdf"), width=9, height=7, units="in", dpi=700)

# plot all simulated trajectories: weighted composition distance only
sim_out |>
  filter(group %in% top_groups$group,
         loss =="weighted_composition_distance") |>
  mutate(year = START_YEAR + year) |>
  ggplot(aes(x=year, colour=loss)) +
  geom_line(aes(y=value, colour=loss, group=interaction(rep,loss)), alpha=0.1) +
  geom_line(data=opt_fit, aes(y=value), linetype="dashed", linewidth=1) +
  facet_wrap(~group, ncol=4, scales="free_y") +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function") +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  clean_theme +
  theme(legend.position="bottom")
ggsave(paste0(here::here("plots"), "/sim_trajectories_weighted-composition.pdf"), width=9, height=7, units="in", dpi=700)



# Combine the above plot with xt and the optimal fits
# ==============================================================================

# data
xt <- read_csv(common$XT_PATH, col_types=cols())

# simulated trajectory summaries
sim_summ <- sim_out |>
  filter(group %in% top_groups$group) |>
  group_by(group, year, loss) |>
  summarise(L = quantile(value, probs=(1-QTILE)/2),
            med = quantile(value, probs=0.5),
            U = quantile(value, probs=1-(1-QTILE)/2),
            .groups = "drop") |>
  mutate(year = START_YEAR + year)

# summarise sensitivity for all 3 loss functions
xt |>
  filter(group %in% top_groups$group) |>
  ggplot(aes(x=year)) +
  geom_point(aes(y=frac), size=1) +
  geom_line(aes(y=frac), linewidth=0.5, linetype="dashed") +
  geom_segment(aes(y=frac-frac_sd*1.96, yend=frac+frac_sd*1.96)) +
  geom_line(data=opt_fit, aes(y=value, colour=loss), linewidth=1) +
  geom_ribbon(data=sim_summ, aes(ymin=L, ymax=U, fill=loss), alpha=0.2) +
  labs(x="Year", y="ST-SC group frequency", colour="Loss function", fill="Loss function",
       title=paste0("Top ",N_GROUPS, " ST-GPSC groups. Ribbons show ", QTILE*100, "% trajectory band.")) +
  scale_colour_brewer(palette="Dark2",
                      labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  scale_fill_brewer(palette="Dark2",
                    labels = c("SSE", "Weighted composition", "Weighted SSE")) +
  facet_wrap(~group, scales="free_y") +
  clean_theme +
  theme(legend.position="bottom")
ggsave(paste0(here::here("plots/"), "/sensitivity_top",N_GROUPS,"_ribbon",QTILE,".pdf"), width=9, height=7, units="in", dpi=700)

# summarise sensitivity for just the chosen loss function
xt |>
  filter(group %in% top_groups$group) |>
  ggplot(aes(x=year)) +
  geom_point(aes(y=frac), size=1) +
  geom_line(aes(y=frac), linewidth=0.5, linetype="dashed") +
  geom_segment(aes(y=frac-frac_sd*1.96, yend=frac+frac_sd*1.96)) +
  geom_line(data=filter(opt_fit, loss==BEST_LOSS), aes(y=value), linewidth=1) +
  geom_ribbon(data=filter(sim_summ, loss==BEST_LOSS), aes(ymin=L, ymax=U), alpha=0.2) +
  labs(x="Year", y="ST-SC group frequency",
       title=paste0("Top ",N_GROUPS, " ST-GPSC groups. Ribbons show ", QTILE*100, "% range."),
       subtitle=paste("Based on the optimal fit from the", BEST_LOSS, "function.")) +
  facet_wrap(~group, scales="free_y") +
  clean_theme +
  theme(legend.position="bottom")
ggsave(paste0(here::here("plots/"), "/sensitivity_top",N_GROUPS,"_ribbon",QTILE,"_",BEST_LOSS,".pdf"), width=9, height=7, units="in", dpi=700)
