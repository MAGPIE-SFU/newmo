# functions for the NFDS-based vaccine design 

compress_vstring = function(vaccine) { paste(vaccine, collapse="-")}

getIPD <- function(simulation, inv, group_lookup, years = 6) {
  cols <- paste0("year.", years)
  missing <- setdiff(cols, names(simulation))
  if (length(missing)) stop("Years not in simulation: ",
                            paste(missing, collapse = ", "))
  
  # simulation needs the serotype, right now the row names are the groups 
  simulation$serotype = group_lookup$serotype[match( rownames(simulation), group_lookup$group)] 
  # mean IPD in the years we are using
  simulation$xmean <- rowMeans(simulation[, cols, drop = FALSE])
  
  dat = simulation %>% select(serotype, xmean) %>% inner_join( inv, by = "serotype") %>% 
    mutate(IPD = xmean*v_child) 

  totalIPD = sum(dat$IPD, na.rm = T) 
  return(list(total= totalIPD, seroIPD = dat %>% select(serotype, IPD)))
  } 

# the years input states the years to average over in computig the IPD score, eg 6:10. 
score_vaccine <- function(modelpars, vseros, inv, group_lookup, years = 6 ) {
  modelpars$vaccine1=vseros
  modelpars$vaccine2 = vseros
  mysim = simulate(samodel, modelpars, nyears=max(years))
  ipd  <- getIPD(mysim, inv, group_lookup, years = years)$total # only cumulative score          
  ipd
}


