# this script pulls the reworked data together and creates the 
# NFDS model structure (a list) for the south africa dataset


G <- as.matrix(read.csv("data/SouthAfrica/Ggrouped.csv", row.names=1))
meta <- read_csv("data/SouthAfrica/meta_clean.csv", show_col_types = FALSE)
locus_by_year <- read.csv("data/SouthAfrica/locus_freqs_deltas.csv", row.names=1, check.names = FALSE)
xt <- read_csv("data/SouthAfrica/xt_grouped.csv", show_col_types = FALSE)

all(rownames(G) == xt$group) # TRUE good 
all(colnames(G) == rownames(locus_by_year)) # TRUE good 

samodel <- list(G =G,
                st = xt$serotype[xt$year == 2008], # serotype for each strain (ie each row of G)
                x0 = xt$count[xt$year == 2008]   ,
                el = as.numeric(locus_by_year$pre2009), # equilibrium frequencies 
                delta = as.numeric(locus_by_year$delta)
) 

# print(samodel)

saveRDS(samodel, file = "data/SouthAfrica/samodel.rds")
