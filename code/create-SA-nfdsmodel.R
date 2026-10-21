# this script pulls the reworked data together and creates the 
# NFDS model structure (a list) for the south africa dataset


Gnfds <- as.matrix(read.csv("data/SouthAfrica/processed-inputs/Gnfds.csv", row.names=1))
meta <- read_csv("data/SouthAfrica/processed-inputs/meta_clean_groups.csv", show_col_types = FALSE)

locus_by_year <- read.csv("data/SouthAfrica/processed-inputs/locus_freqs_deltas.csv", check.names = FALSE)
rownames(locus_by_year) <- locus_by_year$locus      # only if something downstream wants them

xt <- read_csv("data/SouthAfrica/processed-inputs/xt_grouped.csv", show_col_types = FALSE)

x0 <- xt %>%
  filter(year <= 2008) %>%
  group_by(group) %>%
  summarise(frac = mean(frac), .groups = "drop")

x0 <- setNames(x0$frac, x0$group)[rownames(Gnfds)]   # align to G's rows
stopifnot(!anyNA(x0), abs(sum(x0) - 1) < 1e-8)

xt$serotype <- ifelse(xt$group == "rare_NVT", "rare_NVT",
                      sub("_[^_]+$", "", xt$group))

all(rownames(Gnfds) == xt$group) # TRUE good 
all(colnames(Gnfds) == rownames(locus_by_year$locus[locus_by_year$in_model])) # not sure why this matters
# but i guess if you want to make an NFDS model with loci selected on another 
# criterion (z, or logratio, or delta, or whatever) then you can see the relevant info for the 
# relevant columns right in locus_by_year 

samodel <- list(G =Gnfds,
                st = xt$serotype[xt$year == 2008], # serotype for each strain (ie each row of G)
                x0 = x0  , # use mean pre-vaccine and just fit to post-vaccine? 
                el = as.numeric(locus_by_year$pre2009[locus_by_year$in_model]), # equilibrium frequencies 
                logitdiff = as.numeric(locus_by_year$logitdiff)
) 

# print(samodel)

saveRDS(samodel, file = "data/SouthAfrica/samodel.rds")
