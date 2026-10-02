# functional versions of the plots currently in a script (compare_sim_data.R) 


# a couple of helper functions 


# helper: sum frequencies within a grouping variable (e.g. vt or serotype)
agg <- function(df, by, xt) {
  ginfo <- xt |>
    distinct(group, serotype) |>
    mutate(vt = if_else(serotype %in% modelpars$vaccine, "VT", "NVT"))
  
    df |>
    left_join(ginfo, by = "group") |>
    group_by(across(all_of(c(by, intersect(c("year","sim"), names(df)))))) |>
    summarise(frac = sum(frac), .groups = "drop")
}

# plot instructions, this goes in the ggplot call 
compare_layers <- function(d, s, b, colour_var, vax_start, prevax_years) {
  vt_cols <- c(VT = "#D55E00", NVT = "#0072B2")
  list(
    geom_line(data = d, aes(year, frac, colour = .data[[colour_var]]), alpha = 0.3),
    geom_point(data = d, aes(year, frac, colour = .data[[colour_var]]), size = 1.5),
    geom_line(data = s, aes(year, frac, colour = .data[[colour_var]], linetype = sim), linewidth = 1),
    # geom_line(data = s, aes(year, frac, colour = .data[[colour_var]]), linewidth = 1),# pre-list version
    geom_segment(data = b, aes(y = frac, yend = frac, colour = .data[[colour_var]]),
                 x = min(prevax_years), xend = vax_start - 0.5, linetype = "dashed"),
    geom_vline(xintercept = vax_start - 0.5, linetype = "dotted"),
    scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted", "twodash")),
    scale_colour_manual(values = vt_cols), 
    labs(x = NULL, y = "Frequency", colour = NULL, linetype = "Simulation", 
         caption = "Points: data. Lines: simulations (line type = simulation). Dashed: x0 (pre-vaccine mean)."),    
    theme(legend.key.width = unit(1.5, "cm")),
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.6, colour = "black")))
  )
}

make_sim_long <- function(sim, samodel, vax_start) {
  # make sim_long 
  # these have group, fraction and year (and the data one has frac_sd) 
  # model: rows = groups, columns = year.1 ... year.k (make sure each year sums to 1)
  sim_long <- sim |>
    sweep(2, colSums(sim), "/") |>
    as.data.frame() |>
    rownames_to_column("group") |>
    pivot_longer(starts_with("year."), names_to = "k",
                 names_prefix = "year.", values_to = "frac") |>
    mutate(year = vax_start - 1 + as.integer(k)) |>
    select(group, year, frac)
  
  # x0 as frequencies (the pre-vaccine mean the model starts from)
  x0_df <- tibble(group = rownames(samodel$G),
                  frac  = as.numeric(samodel$x0) / sum(samodel$x0))
  bind_rows(x0_df |> mutate(year = vax_start - 0.5), sim_long)
}

# stack a list of simulations -- so we can compare multiple simulations on the same plot 
# sims should be a list , like list(baseline = mysim1, lowerv = mysim2) and so on 
stack_sims <- function(sims, samodel, vax_start) {
  if (is.data.frame(sims)) sims <- list(model = sims)   # allow a single simulation to be passed in, still 
  bind_rows(lapply(sims, make_sim_long, samodel = samodel, vax_start = vax_start), .id = "sim") |>
    mutate(sim = factor(sim, levels = unique(sim)))
}


# plot VTs and NVTs total, model and data together (where data means our synth carriage sample) 
vt_nvt_plot <- function(vax_start=2010, prevax_years = 2000:2009,
                        samodel, xt, sims) { 
  sim_long <- stack_sims(sims, samodel, vax_start)
#  sim_long <- make_sim_long(sims, samodel, vax_start) 
  # group info: serotype, and VT/NVT using the SAME vaccine list as the model
  ginfo <- xt |>
    distinct(group, serotype) |>
    mutate(vt = if_else(serotype %in% modelpars$vaccine, "VT", "NVT"))
  
  data_long <- xt |> select(group, year, frac, frac_sd)
  x0_df <- tibble(group = rownames(samodel$G), frac = as.numeric(samodel$x0) / sum(samodel$x0))
 
  # common layers: data = points (+ faint line), model = thick line,
  # x0 = dashed segment across the pre-vaccine years, dotted line = vaccine start

  # the plot 
  p_vt <- ggplot() +
    compare_layers(agg(data_long, "vt",xt), agg(sim_long, "vt",xt), agg(x0_df, "vt", xt), "vt", vax_start = vax_start,
                   prevax_years = prevax_years) +
    ggtitle("Vaccine-type vs non-vaccine-type frequency")
  return(p_vt)
}


# here's a plot of the top groups 
# n_groups is the number of top ST-SC groups to show 
top_group_plot <- function(vax_start=2010, prevax_years = 2000:2009,
                           n_groups=12, samodel, xt, sims) { 
  # make sim long - this has data, fraction and year (for the groups) 
 # sim_long <- make_sim_long(sims, samodel, vax_start) 
  sim_long <- stack_sims(sims, samodel, vax_start)  # pull out the grouping info (ginfo) 
  ginfo <- xt |>
    distinct(group, serotype) |>
    mutate(vt = if_else(serotype %in% modelpars$vaccine, "VT", "NVT"))
  
  # data long, same thing but with frac_sd as well, from the xt information 
   data_long <- xt |> select(group, year, frac, frac_sd)
   x0_df <- tibble(group = rownames(samodel$G), frac = as.numeric(samodel$x0) / sum(samodel$x0))
   # get the top groups (in the data) 
   top_groups <- data_long |>
    group_by(group) |>
    summarise(m = mean(frac)) |>
    slice_max(m, n = n_groups) |>
    pull(group)

     # look at those groups only 
  d_g <- data_long |> filter(group %in% top_groups) |> left_join(ginfo, by = "group")
  s_g <- sim_long  |> filter(group %in% top_groups) |> left_join(ginfo, by = "group")
  b_g <- x0_df     |> filter(group %in% top_groups) |> left_join(ginfo, by = "group")
  
  # here is the plot , similar to the vt vs nvt plot, but for these grouips. could consider 
  # putting al lthe NVT ones together and all the VT ones together rather than mixing but haven't done that 
  p_groups <- ggplot() +
    geom_errorbar(data = d_g, aes(year, ymin = pmax(frac - frac_sd, 0),
                                  ymax = frac + frac_sd, colour = vt),
                  width = 0, alpha = 0.4) +
    compare_layers(d_g, s_g, b_g, "vt", vax_start = vax_start,
                   prevax_years = prevax_years) +
    facet_wrap(~ group, scales = "free_y") +
    ggtitle(paste("The", n_groups, "most common ST-SC groups"))
  return(p_groups) 
  } 
  

top_sero_plot <- function(vax_start=2010, prevax_years = 2000:2009,
                          n_serotypes=12, samodel, xt, sims ) { 
  # make sim long - this has data, fraction and year (for the groups) 
#   sim_long <- make_sim_long(sims, samodel, vax_start) 
  sim_long <- stack_sims(sims, samodel, vax_start)
  # data long, same thing but with frac_sd as well, from the xt information 
  data_long <- xt |> select(group, year, frac, frac_sd)
  x0_df <- tibble(group = rownames(samodel$G), frac = as.numeric(samodel$x0) / sum(samodel$x0))
  d_st <- agg(data_long, c("serotype", "vt"), xt)
  top_st <- d_st |>
    group_by(serotype) |>
    summarise(m = mean(frac)) |>
    slice_max(m, n = n_serotypes) |>
    pull(serotype)
  
  p_sero <- ggplot() +
    compare_layers(d_st |> filter(serotype %in% top_st),
                   agg(sim_long, c("serotype", "vt"), xt) |> filter(serotype %in% top_st),
                   agg(x0_df, c("serotype", "vt"), xt) |> filter(serotype %in% top_st),
                   "vt", vax_start = vax_start,
                   prevax_years = prevax_years) +
    facet_wrap(~ serotype, scales = "free_y") +
    ggtitle(paste("The", n_serotypes, "most common serotypes"))
  
  p_sero
}
  
  
  
  
  
  
