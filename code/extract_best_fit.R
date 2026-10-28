source("code/simulate_ode.R")
source("code/parameter_generator.R")

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)


rds_params_fit_path <- common$FITTED_PARAMETERS_PATH_SEQ

all_fits <- readRDS(rds_params_fit_path)

fit_summary <- readRDS(file.path(common$FIT_SEQUENTIAL_DIR, "fit_summary_switch.rds"))

idx <- tapply(seq_len(nrow(fit_summary)), fit_summary[, "loss"],
              function(i) i[which.min(fit_summary[i, "fitted"])])
# fit_summary[idx, ]

best <- all_fits[all_fits$fit_id %in% fit_summary[idx, "fit_id"], ]

saveRDS(best, common$BEST_FITS_PATH)

modelpars_list <- setNames(
  Map(function(id, loss, method) {
    parameter_generator(
      how_numerical = loss, optimizer = method, fit_id = id,
      path = common$BEST_FITS_PATH,
      vaccine_type1 = "PCV7", vaccine_type2 = "PCV13",
      how_gene_order = "default",
      use_switch = TRUE
    )
  }, best$fit_id, best$loss, best$method),
  best$fit_id
)

samodel <- readRDS(common$SA_MODEL_PATH)
simulations <- lapply(modelpars_list, function(modelpars) {
    simulate(samodel, modelpars, common$NYEARS)
})

source("code/parameter_fitting_functions.R")

data_sd <- load_fit_data()

# loss_functions[["SSE"]](simulations[[1]], data_sd$data, data_sd$sd)
errors <- t(sapply(simulations, function(sim)
  sapply(loss_functions, function(f) f(sim, data_sd$data, data_sd$sd))
))
rel <- sweep(errors, 2, apply(errors, 2, min), "/")
score <- rowSums(rel)
best_i <- which.min(score)

best_params <- best[best_i, ]

saveRDS(best_params, common$OVERALL_BEST_FITS_PATH)
