# invent vaccines based on a greedy approach, and on an exhaustive search 

# prelims
library(dplyr)
library(patchwork)
library(future.apply)

source("code/simulate_ode.R")

common <- new.env()
sys.source(here::here("code", "common.R"), envir = common)

# Load information
samodel <- readRDS(common$SA_MODEL_PATH)
xt <- read_csv(
  common$XT_PATH,
  show_col_types = FALSE
)
nyears <- common$NYEARS



# ---- method 1: use just the estimated contribution to the IPD population in our inferred dataset 

sero_car <- xt %>%
  filter(year %in% 2000:2008) %>%        # e.g. pre-vaccine only
  group_by(year, serotype) %>%
  summarise(frac = sum(frac), .groups = "drop") %>%   # over groups
  group_by(serotype) %>%
  summarise(frac = mean(frac), .groups = "drop")      # over years


# ---- method 2: greedy 




# ---- method 3: exhaustive search 

