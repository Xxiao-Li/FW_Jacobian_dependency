#============================================================
# 1. Three-species module dependency analysis
#============================================================

rm(list = ls())

library(tidyverse)
library(cowplot)
library(viridis)
library(scales)
library(here)

setwd(here())
set.seed(2026)

n_rep <- 10000
fs  <- 20
fss <- 18

case_levels <- c(
  "(b) Full dependency",
  "(c) c dependency only",
  "(d) B dependency only",
  "(e) e dependency only",
  "(f) No dependency"
)

scenario_labs_short <- c(
  "(b) Full dependency"   = "Full",
  "(c) c dependency only" = "c",
  "(d) B dependency only" = "B",
  "(e) e dependency only" = "e",
  "(f) No dependency"     = "No"
)

# ------------------------------------------------------------------
# Function: compute Re(lambda_max)
# Polynomial: lambda^3 + T1 * lambda + T0 = 0
# ------------------------------------------------------------------

get_re_lammax <- function(t1, t0) {
  coefs <- c(t0, t1, 0, 1)
  roots <- polyroot(coefs)
  max(Re(roots))
}

# ------------------------------------------------------------------
# Function: simulate one dependency scenario
# ------------------------------------------------------------------

simulate_scenario <- function(
    case_name,
    c12_1, c12_2, c13_1, c13_2, c23_1, c23_2,
    B1_1, B1_2, B2_1, B2_2, B3_1, B3_2,
    e2_1, e2_2, e3_1, e3_2
) {
  
  T1 <- (c12_1 * B1_1) * (e2_1 * c12_2 * B2_1) +
    (c13_1 * B1_2) * (e3_1 * c13_2 * B3_1) +
    (c23_1 * B2_2) * (e3_2 * c23_2 * B3_2)
  
  T0 <- (e2_1 * c12_1 * B2_1) *
    (e3_1 * c23_1 * B3_1) *
    (c13_1 * B1_1) -
    (c12_2 * B1_2) *
    (c23_2 * B2_2) *
    (e3_2 * c13_2 * B3_2)
  
  ReLambdaMax <- mapply(get_re_lammax, T1, T0)
  
  tibble(
    T1 = T1,
    T0 = T0,
    ReLambdaMax = ReLambdaMax,
    Case = case_name
  )
}

# ------------------------------------------------------------------
# Scenario (a): analytical landscape
# ------------------------------------------------------------------

t1_seq <- seq(0, 2, length.out = 100)
t0_seq <- seq(-2, 2, length.out = 100)

grid <- expand.grid(
  T1 = t1_seq,
  T0 = t0_seq
) %>%
  as_tibble() %>%
  mutate(
    ReLambdaMax = mapply(get_re_lammax, T1, T0),
    Case = "(a) Analytical"
  )

# ------------------------------------------------------------------
# Scenario (b): full dependency
# ------------------------------------------------------------------

c12 <- runif(n_rep)
c13 <- runif(n_rep)
c23 <- runif(n_rep)

B1 <- runif(n_rep)
B2 <- runif(n_rep)
B3 <- runif(n_rep)

e2 <- runif(n_rep)
e3 <- runif(n_rep)

df_full <- simulate_scenario(
  "(b) Full dependency",
  c12, c12,
  c13, c13,
  c23, c23,
  B1, B1,
  B2, B2,
  B3, B3,
  e2, e2,
  e3, e3
)

# ------------------------------------------------------------------
# Scenario (c): c dependency only
# ------------------------------------------------------------------

c12 <- runif(n_rep)
c13 <- runif(n_rep)
c23 <- runif(n_rep)

df_c <- simulate_scenario(
  "(c) c dependency only",
  c12, c12,
  c13, c13,
  c23, c23,
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep)
)

# ------------------------------------------------------------------
# Scenario (d): B dependency only
# ------------------------------------------------------------------

B1 <- runif(n_rep)
B2 <- runif(n_rep)
B3 <- runif(n_rep)

df_B <- simulate_scenario(
  "(d) B dependency only",
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  B1, B1,
  B2, B2,
  B3, B3,
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep)
)

# ------------------------------------------------------------------
# Scenario (e): e dependency only
# ------------------------------------------------------------------

e2 <- runif(n_rep)
e3 <- runif(n_rep)

df_e <- simulate_scenario(
  "(e) e dependency only",
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  e2, e2,
  e3, e3
)

# ------------------------------------------------------------------
# Scenario (f): no dependency
# ------------------------------------------------------------------

df_no <- simulate_scenario(
  "(f) No dependency",
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep),
  runif(n_rep), runif(n_rep)
)

# ------------------------------------------------------------------
# Combine data
# ------------------------------------------------------------------

df_all <- bind_rows(
  df_full,
  df_c,
  df_B,
  df_e,
  df_no
) %>%
  mutate(
    Case = factor(Case, levels = case_levels)
  )

grid <- grid %>%
  mutate(
    Case = factor(Case, levels = case_levels)
  )

# ------------------------------------------------------------------
# Plot settings
# - analytical panel uses its own wide scale
# - simulated panels use their own range for better contrast
# ------------------------------------------------------------------

fill_limits_analytical <- c(0, 1.4)
fill_limits_sim <- c(0, 0.4)

x_limits_sim <- range(df_all$T1, na.rm = TRUE)
y_limits_sim <- range(df_all$T0, na.rm = TRUE)

base_panel_theme <- theme(
  axis.title = element_text(size = fs),
  axis.text = element_text(size = fss),
  plot.title = element_text(size = fss, face = "bold"),
  plot.margin = margin(3, 3, 3, 3, "pt")
)

# ------------------------------------------------------------------
# Panel (a): analytical landscape
# ------------------------------------------------------------------

p_analytical <- ggplot(grid, aes(x = T1, y = T0, fill = ReLambdaMax)) +
  geom_raster(interpolate = TRUE) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    colour = "grey35"
  ) +
  scale_fill_viridis_c(
    option = "H",
    limits = fill_limits_analytical,
    name = expression(Re(lambda[max]))
  ) +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  labs(
    x = expression(T[1]),
    y = expression(T[0]),
    title = "(a) Analytical"
  ) +
  theme_minimal(base_size = fs) +
  base_panel_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = fss, color = "black"),
    legend.text = element_text(size = fss, color = "black"),
    axis.title = element_text(size = fs, color = "black"),
    axis.text = element_text(size = fss, color = "black"),
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------
# Ratio plot: |T0| / |T1| against Re(lambda_max)
# ------------------------------------------------------------------

grid_plot <- grid %>%
  mutate(
    T0_T1_ratio = ifelse(T1 != 0, abs(T0) / abs(T1), NA_real_),
    T0_sign = ifelse(T0 >= 0, "T0 >= 0", "T0 < 0")
  ) %>%
  filter(is.finite(T0_T1_ratio))

p_ratio <- ggplot(
  grid_plot,
  aes(x = T0_T1_ratio, y = ReLambdaMax, color = T0_sign)
) +
  geom_point(alpha = 0.20, size = 0.4) +
  scale_color_manual(
    values = c("T0 >= 0" = "darkblue", "T0 < 0" = "darkred"),
    labels = c(
      "T0 >= 0" = expression(T[0] >= 0),
      "T0 < 0"  = expression(T[0] < 0)
    ),
    name = expression(sign(T[0]))
  ) +
  scale_x_log10(labels = scales::label_number()) +
  labs(
    x = expression("|" * T[0] * "| / |" * T[1] * "|"),
    y = expression(Re(lambda[max]))
  ) +
  theme_test(base_size = fs) +
  base_panel_theme +
  theme(
    legend.position = c(0.4, 0.96),
    legend.justification = c(1, 1),
    legend.title = element_text(size = fss - 1),
    legend.text = element_text(size = fss - 1),
    axis.title = element_text(size = fs, color = "black"),
    axis.text = element_text(size = fss, color = "black"),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.75),
      colour = NA
    ),
    legend.key = element_rect(fill = NA, colour = NA),
    panel.grid = element_blank()
  )

# ------------------------------------------------------------------
# Function for simulated scenario panels
# ------------------------------------------------------------------

plot_scenario <- function(df, show_legend = FALSE) {
  ggplot(df, aes(x = T1, y = T0, color = ReLambdaMax)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.3,
      colour = "grey35"
    ) +
    geom_point(alpha = 0.45, size = 1) +
    scale_color_viridis_c(
      option = "H",
      limits = fill_limits_sim,
      breaks = c(0, 0.2, 0.4),
      oob = scales::squish,
      name = expression(Re(lambda[max])),
      guide = if (show_legend) {
        guide_colorbar(
          direction = "horizontal",
          barwidth = grid::unit(3.0, "cm"),
          barheight = grid::unit(0.25, "cm"),
          title.position = "top",
          label.position = "bottom"
        )
      } else {
        "none"
      }
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    coord_cartesian(
      xlim = x_limits_sim,
      ylim = y_limits_sim
    ) +
    labs(
      x = expression(T[1]),
      y = expression(T[0]),
      title = unique(as.character(df$Case))
    ) +
    theme_test(base_size = fs) +
    base_panel_theme +
    theme(
      plot.margin = margin(3, 10, 3, 3, "pt"),
      legend.position = if (show_legend) c(0.70, 0.2) else "none",
      legend.justification = c(0.5, 0.5),
      legend.title = element_text(size = fss - 2),
      legend.text = element_text(size = fss - 3),
      axis.title = element_text(size = fs, color = "black"),
      axis.text = element_text(size = fss, color = "black"),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.78),
        colour = NA
      ),
      legend.key = element_rect(fill = NA, colour = NA)
    )
}

# ------------------------------------------------------------------
# Create simulated panels
# ------------------------------------------------------------------

p_full <- plot_scenario(
  df_all %>% filter(Case == "(b) Full dependency"),
  show_legend = FALSE
)

p_c <- plot_scenario(
  df_all %>% filter(Case == "(c) c dependency only"),
  show_legend = FALSE
)

p_B <- plot_scenario(
  df_all %>% filter(Case == "(d) B dependency only"),
  show_legend = FALSE
)

p_e <- plot_scenario(
  df_all %>% filter(Case == "(e) e dependency only"),
  show_legend = FALSE
)

p_no <- plot_scenario(
  df_all %>% filter(Case == "(f) No dependency"),
  show_legend = TRUE
)

# ------------------------------------------------------------------
# Panel (g): scenario comparison
# ------------------------------------------------------------------

scenario_cols <- c(
  "(b) Full dependency"   = "#355C7D",
  "(c) c dependency only" = "#3C8D7B",
  "(d) B dependency only" = "#8C8C39",
  "(e) e dependency only" = "#E39C47",
  "(f) No dependency"     = "#E7A4B0"
)

df_box <- df_all %>%
  mutate(
    Case = factor(Case, levels = rev(case_levels))
  )

p_box <- ggplot(df_box, aes(x = ReLambdaMax, y = Case, fill = Case)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  scale_y_discrete(labels = rev(scenario_labs_short)) +
  coord_cartesian(xlim = c(0, 0.4)) +
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.1)) +
  labs(
    x = expression(Re(lambda[max])),
    y = NULL,
    title = "(g) Scenario comparison"
  ) +
  theme_test(base_size = fs) +
  base_panel_theme +
  theme(
    axis.text.y = element_text(size = fss, color = "black"),
    axis.text.x = element_text(size = fss, color = "black")
  )


p_blank <- ggplot() + theme_void()

row1_inner <- plot_grid(
  p_analytical,
  p_blank,
  p_ratio,
  ncol = 3,
  rel_widths = c(1.1, 0.10, 1.1),
  align = "hv"
)

row1 <- plot_grid(
  p_blank,
  row1_inner,
  p_blank,
  ncol = 3,
  rel_widths = c(0.05, 1, 0.05)
)

row2 <- plot_grid(
  p_full,
  p_c,
  p_B,
  ncol = 3,
  align = "hv"
)

row3 <- plot_grid(
  p_e,
  p_no,
  p_box,
  ncol = 3,
  align = "hv"
)

p_all <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(1.1, 1, 1),
  align = "v"
)

#p_all

ggsave(filename = "figures/3_species_module_dependency.jpeg",
  plot = p_all, width = 14, height = 13, dpi = 600)

ggsave(filename = "figures/3_species_module_dependency.pdf",
  plot = p_all, width = 14, height = 13)


