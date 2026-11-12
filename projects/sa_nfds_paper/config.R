results_dir <- "projects/sa_nfds_paper/results"
META_RAW_PATH <- "raw_data/SouthAfrica/metadata.csv"
INV_CHILDREN_PATH <- "raw_data/Invasiveness/nick-inv-sts-children.xlsx"
INV_ADULTS_PATH   <- "raw_data/Invasiveness/nick-inv-sts-adults.xlsx"


NON_TYPEABLE <- c("untypable", "Swiss_NT")


# Renames: the isolate IS this serotype from here on (applied in clean_metadata)
SEROTYPE_RENAMES <- c(
  "11C"     = "11B",    # Nick has 11B only; group 11B/11C under that name
  "10X"     = "10B",
  "23B1"    = "23B",    # 23B1 is a 23B variant
  "33A/33F" = "33F",
  "15B/15C" = "15B/C",  # Nick's spelling
  "20A"     = "20",     # 20 was split into 20A/20B after Nick's data
  "20B"     = "20",
  "35D"     = "35B"     # 35D isolates were typed as 35B until recently
)

# Invasiveness donors: the serotype stays itself but borrows the donor's estimate
INVASIVENESS_DONORS <- c(
  "11BC" = "11B",
  "6G"   = "6B",
  "24"   = "24F",
  "33D"  = "33A"
)

# Carriage reweighting ----------------------------------------------------------
# sept 22: omega was the weight for carriage compared to IPD; set to 0.27 because that's the
# fraction of the data that's from carriage. Sept 24: replaced by the Dirichlet prior k below
# (better statistical interpretation, works out very similarly).
# Leonie's paper notes NFDS datasets to date only had data from children, so tauA ~ 0.
REWEIGHT_TAU_A <- 0                      # weight for adult IPD, should be v small
REWEIGHT_TAU_U <- 272 / (3865 + 1386)    # prop of the dataset that is unknown age; REDO if new dataset
REWEIGHT_K     <- 100                    # Dirichlet prior strength, in equivalent child carriage swabs

# Group definitions
PREVAX_YEARS  <- 2000:2009
MIN_SERO_FRAC <- 0.001   # NVT serotypes below this share of pre-vaccine synthetic carriage -> rare_NVT
MIN_SG_FRAC   <- 0.005   # serotype-GPSC pairs below this share -> serotype_other

# x(t) uncertainty
# PSEUDOCOUNT was deleted; it was bringing all the _other groups too high and
# pushing things too much towards the mean
XT_LAMBDA <- 0.01   # prior mean smoothed toward uniform; fraction of the prior mean that is uniform
XT_N_REPS <- 1000   # repetitions for the invasiveness uncertainty

GENE_MATRIX_RAW_PATH <- "raw_data/SouthAfrica/pa_matrix_final.csv"

# Locus statistics: pre- and post-vaccine windows
LOCUS_PRE_YEARS  <- 2000:2008
LOCUS_POST_YEARS <- 2012:2014
LOCUS_ALL_YEARS  <- 2000:2014

X0_YEARS <- 2005:2009   # pre-vaccine years averaged for x0

# Model time setup
NYEARS   <- 6
T_SWITCH <- 2 * 12   # 2 years between PCV7 (2009) and PCV13 (2011), 12 time units per year

# Default (CC_MAGIC) parameters: fit starting values and the baseline
DEFAULT_PARS <- list(
  v1      = -log(1 - 0.065),
  sigma   = log(1 + 0.13),
  migr    = 0.00001,
  portion = 0.3,
  whigh   = 0.1363,
  wfrac   = 0.0023 / 0.1363
)

# Sequential fits
FIT_SEEDS         <- c(42, 4224, 1337, 69)
FIT_USE_SWITCH    <- TRUE
FIT_USE_LOG       <- TRUE
FIT_CHAIN_METHODS <- c("nmkb", "lbfgsb", "hjkb")
FIT_N_ROUNDS      <- 3
FIT_START_FROM    <- "screen"
FIT_N_SCREEN      <- 5000
FIT_MAXFEVAL      <- 20000
FIT_REL_TOL       <- NULL     # 1e-6
FIT_PARALLEL      <- TRUE     # parallelise the start screen
FIT_N_WORKERS     <- max(1, parallelly::availableCores() - 1)
# This floor is to avoiod weighting small (rare) groups insanely high
FIT_SD_FLOOR <- 0.004   # minimum frac_sd for the fit weights (median of frac_sd)

# Fit diagnostics
DIAG_N_CORES         <- 1          # > 1 uses parallel::mclapply (Linux/macOS only)
PROFILE_DECADES      <- 1.5        # log-scale parameters: +/- this many decades around the fit
PROFILE_N_GRID       <- 50         # grid points per parameter
PROFILE_YMAX         <- 10         # cap on the loss ratio
PROFILE_RELATIVE_TO  <- "best_seed"  # or "own"
LANDSCAPE_N_SCREEN   <- 1000       # random starts (the default start is added too)
LANDSCAPE_SEED       <- 1337

# Fit figures used to be run_and_plot.R
FIG_SETTINGS <- list(
  vax_start    = 2010,
  prevax_years = 2000:2009,
  n_groups     = 20,
  n_serotypes  = 18,
  vaccine      = nfdsmodelling::VACCINE_TYPES[["PCV13"]]
)

# Sensitivity analysis
SENS_N_GROUPS   <- 16       # number of groups to plot
SENS_QTILE      <- 0.95     # summary quantile for the ribbons
SENS_N_SIMS     <- 100      # parameter draws per loss function
SENS_START_YEAR <- 2007.5   # year prior to start of model evaluation
SENS_BEST_LOSS  <- "SSE"    # loss of the chosen optimal fit, for the last plot
SENS_N_CORES    <- 1        # > 1 uses mclapply (Linux/macOS only)
SENS_SD_PARS    <- data.frame(
  par = c("logv1", "logsigma", "portion"),
  SD  = c( 0.1,     0.25,       0.01),
  # bounds on the sampling scale: portion in [0,1]
  lo  = c(-Inf,    -Inf,        0.01),
  hi  = c( Inf,     Inf,        0.99)
)

# Projected IPD
IPD_N_REPS     <- 1000
IPD_N_CARRIAGE <- 100000
FIRST_SIM_CALENDAR_YEAR <- 2011   # corresponds to year.1 in the simulation
# IPD Plot
IPD_DOT_TOP_N <- 20   # number of serotypes in the IPD dot plot

# Vaccine design
VD_TARGET_SIZE  <- 13         # valency of the invented vaccines
VD_N_CANDIDATES <- 18         # most invasive non-PCV7 serotypes to choose from
VD_SCORE_YEARS  <- 6:10       # simulation years to average the IPD over
VD_N_CORES      <- 1          # > 1 uses mclapply (Linux/macOS only)

# Vaccine heatmap
HEATMAP_VACCINES      <- c("PCV7", "PCV13", "PCV15", "PCV21", "PCV23")   # licensed, in display order
HEATMAP_N_EXHAUSTIVE  <- 6     # exhaustive vaccines to show, lowest IPD first
HEATMAP_METHOD_LABELS <- c(
  greedy       = "Greedy",
  exhaustive   = "Exhaustive",
  ipdfrac      = "todo: ipdfrac",
  mostinvasive = "todo: mostinvasive"
)
