# Notes on model testing 

Here are some notes (from CC) about simulating and testing the NFDS model. 

## 2026-10-02: proposed changes to the model

The `rareNVT` group is a problem. It's defined as those ST-SC combinations that had fewer than 40 observations. That sounds reasonable. I defined them as a group, `rare_NVT`. But consider this: suppose I try the model without NFDS, so setting $\sigma$ to 0. Now, each non-vaccine group should rise basically in proportion to its pre-vaccine frequency among NVTs. That is, if group $i$'s portion of the pre-vaccine NVT population is $f_i = n_i / \sum_{ j \in NVT} n_j$ where $j$ is the group and $n_j$ are the numbers of group $j$ in our synthetic carriage population, then without NFDS, we would expect $i$ to also make a fraction $f_i$ of the post-vaccine population. 
The rare NVT group had about 7.5% of the total carriage population pre-vaccine. The NVT strains together were about 33% of the pre-vaccine population. So without NFDS, we would expect the 7.5% to become 7.5/33 or about 23%. That is what happens with $\sigma=0$ but of course, *by definition* the rare NVT group can't rise in the data because the strain is conditioned to stay rare. Now, of course, we will have NFDS in our model, so the no-NFDS thought experiment isn't the same as the simulations we will have. But even with NFDS I was seeing large rises in rareNVTs which will never be possible in the data, and I was not seeing rises in some key NVTs. 


*Solution*" This is a setup issue - I will re-create the whole model and $x(t)$ and so on with a different definition of rare NVTs, using only pre-vaccine data to decide. This will mean we'll have more than the current 87 groups, but I think it's worth it. 


The NVT strains that do rise (in the data) may rise only very slowly in the model if migration is very low and they start at low levels in the initial conditions. We might want to start strains that have near-0 frequencies pre-vaccine at a somewhat higher (than 0) starting point, and/or change the migration dynamics. 

Currently, if the migration is too high, it pulls the system towards a more uniform distribution of strains. Should we have vaccine strains not migrating in as much (since the population is immune to them, they may not establish, as well as failing to reproduce)? Should we have rare strains not migrating in as much as more common strains? (such that the equilibrium dynamics of the model are really the pre-vaccine frequencies? with constant migration they are not quite. But previous models used constant migration). 

Next steps after these changes: 

- make a quantitative fit assessment (e.g. a likelihood or loss function) so that we can fit the model.  Our loss function should probably reflect the fact that our data are proportions, and should probably not end up hugely over-focused on only the most common few strains. 
- consider whether we really need 4 different parameters that are all about how much NFDS there is. More parameters: harder fitting. 
- fit the model parameters
- is there evidence that the fits under NFDS are better than without NFDS? (if there were one NFDS parameter, is there statistical evidence that it differs from 0?) 

Ideally we might also propagate uncertainty better into $x(t)$ than we do now - we have something derived from the numbers of samples but not from the uncertainty in invasiveness, which is considerable. 


## 2026-09-29: plotters and baseline parameters 

I decided that the ode solver should give us not the values of the group frequencies
(or counts, if we are using a population of $N$) but instead, the average for the given
year. Since our synthetic carriage population $x(t)$ is the average for each year, this
seems fair to the model. I changed this in the function `simulate` in `simulate_ode.R`. 


I made three plotters (with help from AI, but I checked every line and recreated them as functions) and their associated helper functions. 
These are in `plot_functions.R` and the script I'm using them in is `cc-notes-tmp.R`, where I am exploring the model simulations
and how they depend on the parameters.
The three plotters are: 

- `vt_nvt_plot(vax_start=2010, prevax_years = 2000:2009,
              samodel, xt, sims)` : takes in the model, the synthetic carriage population ("data") $x(t)$ and sims, a list of simulations with names.
  This makes a plot showing the VT and NVT totals in the simulation, compared to the data. These are, naturally, most sensitive to the model parameter $v$, the
  strength of the vaccine effect. However, changing the NFDS parameters does change these eventually. More loci under NFDS, a stronger NFDS effect and a higher parameter `wlow` all mean
  that VT strains aren't reduced as fast as they are in the data and NVT strains don't rise as fast (at least, in my preliminary exploration).

- `top_group_plot(vax_start=2010, prevax_years = 2000:2009,
 n_groups=18, samodel, xt, sims)`: same idea as the VT/NVT plot, but for the top ST-SC groups. These are the top groups in the data, according to their mean fraction over the years.

- `top_sero_plot <- function(vax_start=2010, prevax_years = 2000:2009,
                          n_serotypes=12, samodel, xt, sims )`: as for top ST-SC groups, but now grouped by serotype (so, all groups with ST 19A together). 

Note that in the first version of my script, PCV13 was mis-specified and should contain 19A. Have updated the model parameters (currently just in my notes script). 

Did some playing around with parameters. We might want to avoid a group like rareNVT that has lots of loci -- it can rise to fill any gap when strains are removed by the vaccine. One thing
I have noticed is that we tend to have too much rareNVT and too little of some specific NVT serotypes (at baseline, 15A, 23A and 35B). 

Here are the baseline parameters, taken from previous NFDS models (eg my MA parameters). Rates are per month, and are converted in the simulation code. I have for now been using a very low migration rate.
But note that the previous value would have been in a population of size 100K but here we have frequencies so the sum is 1. anyway, too high a migration rate means that the VTs don't decline nearly as much
as the data, presumably because they are introduced. Perhaps we should model migration but not the subsequent success of vaccine strains (since our population is immune to them). 

```
 modelpars <-list()
 modelpars$v <- 0.065 # there is a reference to 0.1 in my notes, too..
 modelpars$sigma <- 0.13 # from previous NFDS models
 modelpars$migr <- 0.00001 # migration parameter
 modelpars$portion <- 0.3 # fraction under NFDS
 modelpars$vaccine = c("1", "3", "4", "5", "6A","6B","7F","9V","14","18C","19A", "19F","23F") # PCV13 i think
 modelpars$whigh = 0.1363 # from MA model and previous paper
 modelpars$wlow = 0.0023 # from MA model too; could be 0 (Leonie)

 modelpars$K=sum(samodel$x0) # that's 1.. (often we'd use a carrying capacicty but i don't think it matters and x(t) is frequencies.. )
```

Note that there are *four* distinct parameters that all essentially increase or decrease the amount of NFDS that is going on. `sigma` is the strength of NFDS in the ODE (in the term $(1-\sigma)^\pi$). Increasing it makes the NDFS push the strains up or down more. The `portion` parameter is the fraction of the loci that have a high NFDS weight. The two $w$s `whigh` and `wlow` are the weights given to the loci, in the sum that creates $\pi$.  It might be worth reducing these to
only two parameters: a fraction of genes and a strength. We could keep `wlow` at 0 (I think that is what Leonie did) or some low number, or even have a third parameter for the relative high vs low numbers or something like that. 
Fitting fewer parameters will be easier than fitting more. 






