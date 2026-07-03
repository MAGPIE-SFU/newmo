# Pneumococcus NFDS modelling

This is a private repo for MAGPIE's pneumococcus modelling. 

## Data 
The data/ folder has USA-Mass, SouthAfrica and Invasiveness folders. The Invasiveness folder contains the serotype-specific invasiveness values from Nick Croucher's 2022 paper (Disease progression rates... ). The other two contain the data for Massachusetts and South Africa. 

## Code 
The code/ folder contains scripts that group and reweight the original data to create "the NFDS model" (for South Africa). 

## An NFDS Model 

We refer to the information needed to define the NFDS dynamics (for a location) as "the NFDS model" or "the model" (for that location).

For now, we refer to the gpsc-st groups as "strains". 

The information needed to simulate an NFDS model forward in time is
 
1. G $G_{jl}$ is the fraction of isolates in strain $j$ that have the gene at locus $l$. The columns in order of most stable loci first. 
1. The serotype for each strain (this is also in the strain name, given in the rownames of G).
1. el (the equilibrium frequencies).
1. x0 (the starting values of the populations for each strain).

For the South Africa data (published in Bellman et al), this is in the form of a list in a file called samodel.rds, in the data/ folder.
  
There is a file called simulate_ode.R which does what it says: simulates the ode version of the NFDS dynamics.
