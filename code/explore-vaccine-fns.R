# functions for the NFDS-based vaccine design 


getIPD <- function(simulation, inv,group_lookup) {
  # simulation needs the serotype, right now the row names are the groups 
  simulation$serotype = group_lookup$serotype[match( rownames(simulation), group_lookup$group)] 
  dat = simulation %>% select(serotype, year.6) %>% inner_join( inv, by = "serotype") %>% 
    mutate(IPD = year.6*v_child) 

  totalIPD = sum(dat$IPD, na.rm = T) 
  return(list(total= totalIPD, seroIPD = dat %>% select(serotype, IPD)))
  } 

score_vaccine <- function(modelpars, vseros, inv, group_lookup) {
  modelpars$vaccine1=vseros
  modelpars$vaccine2 = vseros
  mysim = simulate(samodel, modelpars, nyears=6)
  ipd  <- getIPD(mysim, inv, group_lookup)$total # only cumulative score          
  ipd
}


