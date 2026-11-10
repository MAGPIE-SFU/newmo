results_dir <- "projects/sa_nfds_paper/results"
META_RAW_PATH <- "raw_data/SouthAfrica/metadata.csv"
INV_CHILDREN_PATH <- "raw_data/Invasiveness/nick-inv-sts-children.xlsx"
INV_ADULTS_PATH   <- "raw_data/Invasiveness/nick-inv-sts-adults.xlsx"


NON_TYPEABLE <- c("untypable", "Swiss_NT")


# Hard-coded substitutions; these take priority over exact matches and serogroup averages
SEROTYPE_SUBSTITUTIONS <- c(
  "15B/15C" = "15B/C",   # label mismatch only — exact donor exists
  "11BC"    = "11B",
  "35D"     = "35B",     # wciG-deficient 35B variant
  "6G"      = "6B",
  "24"      = "24F",
  "33D"     = "33A"
)
