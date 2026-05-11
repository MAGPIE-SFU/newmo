library(deSolve)

# minimal pneumo model for massachusetts following the earlier matlab version

# read the data 
G <- read.csv( "G_massachussetts.csv", row.names = 1)
metadata <- read.csv("metadata_massachussetts.csv", row.names = 1 ) 
init <- read.csv("initialconditions_mass.csv",row.names = 1)

sttime <- read.csv("STcounts_over_time.csv", header = TRUE)
# these are the counts for the serotypes at the 3 times points, 0, 36 months, 72 months 
# from the data 

# some information for the model 
datast <- as.matrix(sttime[, 2:4]) 

modelpars <-list()
modelpars$v <- 0.061 # there is a reference to 0.1 in my notes, too.. 
modelpars$sigma <- 0.18 
modelpars$migr <- 0.001 # i think actually i had m=0 in my ode version before.. ?? 

locweights <- read.csv("locusweights_mass.csv", row.names = 1) # 
# weighteddev = as.matrix(read.csv("weighteddev_mass.csv",row.names = 1)) #  not needed ; see below

# ---- sanity check for frequencies of the loci at the start 
# ifreq = read_tsv("~/Dropbox/Dynamics/NickCroucher/OptimVaccine/FullStrains/RoundOne/mass.locusFreq.csv",
#                col_names = c("locus","freq"))
# locweights <- locweights %>%  left_join(ifreq, by = "locus")
# are these the locus frequencies i get from the initial conditions? 
thesefs <- (t(as.numeric(init$init)) %*% as.matrix(G)) / sum(init$init) 
# plot(locweights$freq, thesefs)
# they are exactly the same. ok, i do not need to load the frequencies from the previous files. 
# 
locweights$freq <- t(thesefs) # these are the initial locus frequencies 
# and they were also used as the equilibrium locus frequencies in the math, I believe,
# as long as the values in 'init' were the initial levels of the strains 

pcv7 <- c("4","6A","6B","9V","14","18C","19F","23F") # the serotypes in the PCV7 vaccines 
seronames <- sttime$names 
serovt <- rep(0, length(seronames)) 
serovt[seronames %in% pcv7] <- 1  # not sure if i use this now? 

VT <- rep(0, nrow(G)) 
VT[ metadata$serotype %in% pcv7] <- 1 # this is whether each strain is in the vaccine 
VT <- VT*modelpars$v   # now it is the v parameter for each strain 


# now set up the model ODE 

rr  <- -log(1 - VT)
rho <- log(1 + modelpars$sigma)
migr <- modelpars$migr # 0.001: a small rate compared to VT (0.06 for vaccine strains) 

K      <- 100000 # populatoin size (carrying capacity) 
LL     <- nrow(G)
tmax   <- 6 # years 
tspan  <- seq(0, tmax * 12, by = 12) # time steps are monthly 


n_loci <- nrow(locweights)

modifiedGT  <- t(G) * locweights$weight # this is 1090 x 603
GtimesModG  <- as.matrix(G) %*% modifiedGT

# we used something called 'weighteddev' -- these were (from getNickData.m):
# result.weighteddev=result.G*(locusweights.*locusfreq)';  
# 
# the reason to structure the computation of pi_i using these is that it minimizes 
# the number and size of matrix multiplications.
# instead of computing f_l using t(x) %*% G/sum(x) 
# and then computing e_l - f_l and then computing G %*% (weights * (e_l - f_l)), 
# we compute G %*% weights * e_l once only, and compress the frequency part 

weighteddev = as.matrix(G) %*% (as.numeric(locweights$weight * locweights$freq))
# since we are re-doing this, i might name these more sensible things this time, i suppose.. 
# note in my version of the model, I was not changing or fitting the weight (locus weight e_l)
# in our upcoming work we will want to be able to do that. 
# but we should still avoid doing the larger matrix multiplications at each dt or each iteration
# and instead use the (new) weights to define matrices for the model just once before the simulation starts 

# anyway here is the ODE that propagates the population forward 
# i guess we need to link parms to modelpars .. but anyway .. 
diffxt <- function(t, x, parms) {
    totalx <- sum(x)
    dx <- x * (
    log(K / totalx) -
      rr +
      rho * (weighteddev - GtimesModG %*% (x / totalx))
  ) + migr
  list(as.numeric(dx))
}

# and here is how to run it 
out <- ode(
  y = as.numeric(init$init),
  times = tspan,
  func = diffxt,
  parms = NULL,
  method = "lsoda",
  atol = 1e-8,
  rtol = 1e-5
)
# the output is 6 rows of 604 numbers
# the first column are the times 
# each row is X(t) for that time
# that is a very counterintuitive format 
vals <- t(out[, -1, drop = FALSE])
colnames(vals) <- out[, 1]
outdf <- as.data.frame(vals); rm(vals)

# in outdf, each column corresponds to a time, and there are 603 rows, one for each strain. 
# from this it should be straightforward to compute the ST changes over time and compare to the data 
head(outdf)

