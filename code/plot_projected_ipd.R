# ============================================================
# Predicted South African IPD:
# Serotype heatmap and uncertainty dot plot
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(forcats)
library(scales)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

INPUT <- "data/SouthAfrica/processed-inputs/simulated_ipd_serotype_uncertainty.csv"
OUT_DIR <- "plots/predictedIPD"

TOP_N <- 20 # Number of serotypes in dot plot

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load predicted IPD composition
# ------------------------------------------------------------

ipd_serotype_summary <- read_csv(INPUT, show_col_types = FALSE)

required <- c(
  "year", "serotype", "ipd_mean",
  "frac", "frac_lower", "frac_upper"
)

missing_cols <- setdiff(required, names(ipd_serotype_summary))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

ipd_serotype_summary <- ipd_serotype_summary %>%
  mutate(
    serotype = as.character(serotype),
    year = as.integer(year)
  ) %>%
  filter(
    !is.na(serotype),
    !is.na(year),
    is.finite(ipd_mean),
    !serotype %in% c("Swiss_NT", "untypable")
  )

# ------------------------------------------------------------
# 2. Order serotypes by overall predicted IPD burden
# ------------------------------------------------------------

sero_order <- ipd_serotype_summary %>%
  group_by(serotype) %>%
  summarise(
    total = sum(ipd_mean),
    .groups = "drop"
  ) %>%
  arrange(desc(total)) %>%
  pull(serotype)

ipd_serotype_summary <- ipd_serotype_summary %>%
  mutate(
    serotype = factor(
      serotype,
      levels = rev(sero_order)
    )
  )

# ============================================================
# Panel A: Total predicted IPD with serotype contributions
# ============================================================



# Identify serotypes contributing the most IPD overall
top_seros <- ipd_serotype_summary %>%
  group_by(serotype) %>%
  summarise(
    total_ipd = sum(ipd_mean),
    .groups = "drop"
  ) %>%
  pull(serotype)

# Prepare stacked IPD data
ipd_stacked <- ipd_serotype_summary %>%
  mutate(
    serotype_plot = if_else(
      serotype %in% top_seros,
      as.character(serotype),
      "Other"
    )
  ) %>%
  group_by(year, serotype_plot) %>%
  summarise(
    ipd_mean = sum(ipd_mean),
    .groups = "drop"
  )

# Order serotypes by total contribution
sero_order <- ipd_stacked %>%
  group_by(serotype_plot) %>%
  summarise(
    total = sum(ipd_mean),
    .groups = "drop"
  ) %>%
  arrange(desc(total)) %>%
  pull(serotype_plot)

ipd_stacked <- ipd_stacked %>%
  mutate(
    serotype_plot = factor(
      serotype_plot,
      levels = sero_order
    )
  )

# Calculate the total predicted IPD for each year
total_ipd <- ipd_stacked %>%
  group_by(year) %>%
  summarise(
    total = sum(ipd_mean),
    .groups = "drop"
  )

# ============================================================
# Plot
# ============================================================

p_total_ipd <- ggplot() +
  geom_area(
    data = ipd_stacked,
    aes(
      x = year,
      y = ipd_mean,
      fill = serotype_plot
    ),
    position = "stack",
    alpha = 0.85
  ) +
  geom_line(
    data = total_ipd,
    aes(x = year, y = total),
    colour = "black",
    linewidth = 1.1
  ) +
  geom_point(
    data = total_ipd,
    aes(x = year, y = total),
    colour = "black",
    size = 2
  ) +
  scale_fill_viridis_d(
    option = "mako",
    begin = 0.15,
    end = 0.9
  ) +
  scale_x_continuous(
    breaks = sort(unique(ipd_stacked$year))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Simulation year",
    y = "Expected IPD cases per 100,000 carriage events",
    fill = "Serotype",
    title = "Predicted invasive pneumococcal disease",
    subtitle = "Total IPD and serotype-specific contributions"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p_total_ipd)


# ============================================================
# PANEL B: Dot plot with 95% uncertainty intervals
# ============================================================

# Identify top serotypes by total expected IPD
sero_order <- ipd_serotype_summary %>%
  group_by(serotype) %>%
  summarise(
    total_ipd = sum(mean_ipd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_ipd)) %>%
  pull(serotype)

top_seros <- head(sero_order, TOP_N)

# Prepare data for Panel B
dot_data <- ipd_serotype_summary %>%
  filter(serotype %in% top_seros) %>%
  mutate(
    serotype = factor(serotype, levels = rev(top_seros)),
    year = factor(year, levels = sort(unique(year)))
  )


# Intervals represent uncertainty in invasiveness,
# not Poisson predictive intervals.

# All years in one panel, coloured and dodged by year


dot_data <- dot_data %>%
  mutate(
    year = factor(year, levels = sort(unique(year)))
  )

dodge <- position_dodge(
  width = 0.75,
  reverse = TRUE
)

p_dot <- ggplot(
  dot_data,
  aes(
    x = 100 * frac,
    y = serotype,
    colour = year,
    group = year
  )
) +
  geom_errorbar(
    aes(
      xmin = 100 * frac_lower,
      xmax = 100 * frac_upper
    ),
    orientation = "y",
    position = dodge,
    width = 0.4,
    linewidth = 0.6,
    alpha = 0.35
  ) +
  geom_point(
    position = dodge,
    size = 2.2,
    alpha = 0.9
  ) +
  scale_colour_viridis_d(
    option = "mako",
    begin = 0.8,
    end = 0.1,
    direction = 1
  ) +
  scale_x_continuous(
    labels = scales::label_number(suffix = "%"),
    limits = c(0, NA)
  ) +
  labs(
    x = "Proportion of total predicted IPD (%)",
    y = "Serotype",
    colour = "Simulation year",
    title = "Predicted serotype contribution to IPD",
    subtitle = paste0(
      "Top ", length(top_seros),
      " serotypes by total expected IPD"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p_dot)


# ============================================================
# 3. Save figures
# ============================================================

ipd_height <- max(
  6,
  2 + 0.18 * length(sero_order)
)

dot_height <- max(
  6,
  3 + 0.31 * length(top_seros)
)

ggsave(
  file.path(OUT_DIR, "ipd_serotype_dotplot.pdf"),
  p_dot,
  width = 11,
  height = dot_height,
  limitsize = FALSE
)

ggsave(
  file.path(OUT_DIR, "ipd_serotype_total.pdf"),
  p_total_ipd,
  width = 11,
  height = ipd_height,
  limitsize = FALSE
)

message("Plots saved to: ", OUT_DIR)
