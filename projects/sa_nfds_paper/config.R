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
