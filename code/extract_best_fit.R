source("code/simulate_ode.R")
source("code/parameter_generator.R")
source(here::here("code", "parameter_fitting_functions.R"))
library(dplyr)

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

# ---- Find all seed runs -------------------------------------------------------
seed_dirs <- list.dirs(common$FIT_SEQUENTIAL_DIR, recursive = FALSE)
seed_dirs <- seed_dirs[grepl("^seed_", basename(seed_dirs))]
stopifnot("No seed_* folders found" = length(seed_dirs) > 0)

read_seed <- function(dir, prefix) {
  f <- file.path(dir, paste0(prefix, basename(dir), ".rds"))
  if (!file.exists(f)) {
    message("Skipping ", basename(dir), ": ", basename(f), " not found (run unfinished?)")
    return(NULL)
  }
  readRDS(f)
}

summary_all <- bind_rows(lapply(seed_dirs, read_seed, prefix = "summary_with_start_"))
params_all  <- bind_rows(lapply(seed_dirs, read_seed, prefix = "parameters_with_start_"))

# parameter_generator reads fits from a file; point each seed at its own
param_paths <- setNames(
  file.path(seed_dirs, paste0("parameters_with_start_", basename(seed_dirs), ".rds")),
  sub("^seed_", "", basename(seed_dirs))
)

# ---- Best fit per seed x loss (excluding the round_00 starting points) --------
best <- summary_all |>
  filter(method != "start", is.finite(fitted)) |>
  group_by(seed, loss) |>
  slice_min(fitted, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(seed, fit_id, loss, method, start_id, fitted) |>
  left_join(params_all |> select(-loss, -method, -start_id),
            by = c("seed", "fit_id")) |>
  mutate(key = paste0("seed", seed, "_", fit_id))   # fit_id repeats across seeds

saveRDS(best, common$ALL_BEST_FITS_PATH)

# ---- Simulate every candidate --------------------------------------------------
modelpars_list <- setNames(lapply(seq_len(nrow(best)), function(i) {
  parameter_generator(
    how_numerical = best$loss[i], optimizer = best$method[i],
    fit_id = best$fit_id[i], path = param_paths[[as.character(best$seed[i])]],
    vaccine_type1 = "PCV7", vaccine_type2 = "PCV13",
    how_gene_order = "default", use_switch = TRUE
  )
}), best$key)

samodel <- readRDS(common$SA_MODEL_PATH)
simulations <- lapply(modelpars_list, simulate,
                      samodel = samodel, nyears = common$NYEARS)

# ---- Cross-evaluate every candidate on every loss -----------------------------
data_sd <- load_fit_data()
errors <- t(sapply(simulations, function(sim)
  sapply(loss_functions, function(f) f(sim, data_sd$data, data_sd$sd))
))
rel <- sweep(errors, 2, apply(errors, 2, min), "/")   # relative to best across ALL seeds
best$score <- rowSums(rel)[best$key]

# ---- Report --------------------------------------------------------------------
cat("\nScore = sum over loss functions of (loss / best loss across all fits);",
    "lower is better, minimum possible =", ncol(errors), "\n")

for (s in sort(unique(best$seed))) {
  tab <- best |>
    filter(seed == s) |>
    arrange(score) |>
    select(loss, method, start_id, fitted, score)
  cat(sprintf("\nSeed %s  ->  best: %s (score %.3f)\n", s, tab$loss[1], tab$score[1]))
  print(as.data.frame(tab), row.names = FALSE, digits = 4)
}

best_params <- best |> slice_min(score, n = 1, with_ties = FALSE)
cat(sprintf("\nOverall best: seed %s, fitted with %s (%s, %s), score %.3f\n",
            best_params$seed, best_params$loss, best_params$method,
            best_params$start_id, best_params$score))

cat("\nWins per loss function across seeds:\n")
print(best |> group_by(seed) |> slice_min(score, n = 1, with_ties = FALSE) |>
        ungroup() |> count(loss, name = "wins"))

saveRDS(best_params, common$OVERALL_BEST_FITS_PATH)