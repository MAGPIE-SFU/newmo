library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# this script assumes that you have samodel which has its x0, you have the simulation (here called mysim), etc

# LB Code to get the values read in to make the plots...
# mysim <- sim
# samodel <- get_samodel(sa_model_path)

# modelpars <- as.list(unlist(DEFAULT_PARS)[PAR_NAMES])

# pars <- DEFAULT_PARS
# pars[names(modelpars)] <- modelpars

#   # wfrac = wlow / whigh, so wlow <= whigh always holds (no clamp, no dead zone).
# if ("wfrac" %in% names(pars)) {
#     pars$wlow  <- pars$whigh * pars$wfrac
#     pars$wfrac <- NULL
#   }

# pars <- c(pars, FIXED_PARS, list(K = sum(samodel$x0)))


# ---- plot settings - only for vis, this doesn't affect the simulation  -------------------------------------------------------------
vax_start    <- 2010          # sim year.1 = 2009, ..., year.6 = 2014
prevax_years <- 2000:2009     # years averaged to make x0
n_groups     <- 12            # how many common ST-SC groups to show
n_serotypes  <- 12            # how many common serotypes to show

# ---- 1. put model and data in the same long format -------------------------

# model: rows = groups, columns = year.1 ... year.k (make sure each year sums to 1)
sim_long <- mysim |>
    sweep(2, colSums(mysim), "/") |>
    as.data.frame() |>
    rownames_to_column("group") |>
    pivot_longer(starts_with("year."), names_to = "k",
                 names_prefix = "year.", values_to = "frac") |>
    mutate(year = vax_start - 1 + as.integer(k)) |>
    select(group, year, frac)

# x0 as frequencies (the pre-vaccine mean the model starts from)
x0_df <- tibble(group = rownames(samodel$G),
                frac  = as.numeric(samodel$x0) / sum(samodel$x0))
sim_long <- bind_rows(x0_df |> mutate(year = vax_start - 0.5), sim_long)
# group info: serotype, and VT/NVT using the SAME vaccine list as the model
ginfo <- xt |>
    distinct(group, serotype) |>
    mutate(vt = if_else(serotype %in% modelpars$vaccine, "VT", "NVT"))

data_long <- xt |> select(group, year, frac, frac_sd)

# helper: sum frequencies within a grouping variable (e.g. vt or serotype)
agg <- function(df, by) {
    df |>
        left_join(ginfo, by = "group") |>
        group_by(across(all_of(c(by, intersect("year", names(df)))))) |>
        summarise(frac = sum(frac), .groups = "drop")
}

vt_cols <- c(VT = "#D55E00", NVT = "#0072B2")

# common layers: data = points (+ faint line), model = thick line,
# x0 = dashed segment across the pre-vaccine years, dotted line = vaccine start
compare_layers <- function(d, s, b, colour_var) {
    list(
        geom_line(data = d, aes(year, frac, colour = .data[[colour_var]]), alpha = 0.3),
        geom_point(data = d, aes(year, frac, colour = .data[[colour_var]]), size = 1.5),
        geom_line(data = s, aes(year, frac, colour = .data[[colour_var]]), linewidth = 1),
        geom_segment(data = b, aes(y = frac, yend = frac, colour = .data[[colour_var]]),
                     x = min(prevax_years), xend = vax_start - 0.5, linetype = "dashed"),
        geom_vline(xintercept = vax_start - 0.5, linetype = "dotted"),
        scale_colour_manual(values = vt_cols),
        labs(x = NULL, y = "Frequency", colour = NULL,
             caption = "Points: data. Solid line: model (annual average). Dashed: x0 (pre-vaccine mean)."),
        theme_bw()
    )
}

# ---- 2. VT vs NVT totals ---------------------------------------------------
p_vt <- ggplot() +
    compare_layers(agg(data_long, "vt"), agg(sim_long, "vt"), agg(x0_df, "vt"), "vt") +
    ggtitle("Vaccine-type vs non-vaccine-type frequency")
p_vt


# ---- 3. common ST-SC groups --------------------------------------------------
top_groups <- data_long |>
    group_by(group) |>
    summarise(m = mean(frac)) |>
    slice_max(m, n = n_groups) |>
    pull(group)

d_g <- data_long |> filter(group %in% top_groups) |> left_join(ginfo, by = "group")
s_g <- sim_long  |> filter(group %in% top_groups) |> left_join(ginfo, by = "group")
b_g <- x0_df     |> filter(group %in% top_groups) |> left_join(ginfo, by = "group")

p_groups <- ggplot() +
    geom_errorbar(data = d_g, aes(year, ymin = pmax(frac - frac_sd, 0),
                                  ymax = frac + frac_sd, colour = vt),
                  width = 0, alpha = 0.4) +
    compare_layers(d_g, s_g, b_g, "vt") +
    facet_wrap(~ group, scales = "free_y") +
    ggtitle(paste("The", n_groups, "most common ST-SC groups"))
p_groups



# ---- 4. common serotypes ---------------------------------------------------
d_st <- agg(data_long, c("serotype", "vt"))
top_st <- d_st |>
    group_by(serotype) |>
    summarise(m = mean(frac)) |>
    slice_max(m, n = n_serotypes) |>
    pull(serotype)

p_sero <- ggplot() +
    compare_layers(d_st |> filter(serotype %in% top_st),
                   agg(sim_long, c("serotype", "vt")) |> filter(serotype %in% top_st),
                   agg(x0_df, c("serotype", "vt")) |> filter(serotype %in% top_st),
                   "vt") +
    facet_wrap(~ serotype, scales = "free_y") +
    ggtitle(paste("The", n_serotypes, "most common serotypes"))

p_sero

# ---- 5. all groups at once: model vs data, one panel per year ---------------

eps <- 1e-4   # so zeros show on the log scale
scat <- inner_join(data_long, sim_long, by = c("group", "year"),
                   suffix = c("_data", "_model")) |>
    left_join(ginfo, by = "group")

p_scatter <- ggplot(scat, aes(frac_data + eps, frac_model + eps, colour = vt)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = vt_cols) +
    facet_wrap(~ year) +
    labs(x = "Data frequency", y = "Model frequency", colour = NULL,
         title = "All ST-SC groups: model vs data") +
    theme_bw()

p_scatter

 # print(p_vt); print(p_groups); print(p_sero); print(p_scatter)
