# =============================================================================
# Three-species module dependency analysis
# =============================================================================

rm(list = ls())

library(tidyverse)
library(cowplot)
library(viridis)
library(scales)
library(here)
library(emmeans)
library(multcomp)
library(multcompView)

set.seed(2026)

n_rep <- 10000
fs <- 20
fss <- 18

case_levels <- c(
  "(b) Full dependency",
  "(c) Feeding-coefficient\ndependency",
  "(d) Shared-biomass\ndependency",
  "(e) Shared-efficiency\ndependency",
  "(f) No dependency"
)

scenario_labs_short <- c(
  "(b) Full dependency" = "Full",
  "(c) Feeding-coefficient\ndependency" = "Feeding coefficient",
  "(d) Shared-biomass\ndependency" = "Shared biomass",
  "(e) Shared-efficiency\ndependency" = "Shared efficiency",
  "(f) No dependency" = "No"
)

scenario_cols <- c(
  "(b) Full dependency" = "#1B9E77",
  "(c) Feeding-coefficient\ndependency" = "#7570B3",
  "(d) Shared-biomass\ndependency" = "#E6AB02",
  "(e) Shared-efficiency\ndependency" = "#D95F02",
  "(f) No dependency" = "#666666"
)

get_re_lmax <- function(T1, T0) {
  max(Re(polyroot(c(T0, T1, 0, 1))))
}

simulate_scenario <- function(case, share_c, share_B, share_e) {
  c12 <- runif(n_rep)
  c13 <- runif(n_rep)
  c23 <- runif(n_rep)
  B1 <- runif(n_rep)
  B2 <- runif(n_rep)
  B3 <- runif(n_rep)
  e2 <- runif(n_rep)
  e3 <- runif(n_rep)
  
  c12a <- c12
  c13a <- c13
  c23a <- c23
  B1a <- B1
  B2a <- B2
  B3a <- B3
  e2a <- e2
  e3a <- e3
  
  c12b <- if (share_c) c12a else runif(n_rep)
  c13b <- if (share_c) c13a else runif(n_rep)
  c23b <- if (share_c) c23a else runif(n_rep)
  
  B1b <- if (share_B) B1a else runif(n_rep)
  B2b <- if (share_B) B2a else runif(n_rep)
  B3b <- if (share_B) B3a else runif(n_rep)
  
  e2b <- if (share_e) e2a else runif(n_rep)
  e3b <- if (share_e) e3a else runif(n_rep)
  
  T1 <- (c12a * B1a) * (e2a * c12b * B2a) +
    (c13a * B1b) * (e3a * c13b * B3a) +
    (c23a * B2b) * (e3b * c23b * B3b)
  
  T0 <- (e2a * c12a * B2a) *
    (e3a * c23a * B3a) *
    (c13a * B1a) -
    (c12b * B1b) *
    (c23b * B2b) *
    (e3b * c13b * B3b)
  
  tibble(
    T1,
    T0,
    ReLambdaMax = map2_dbl(T1, T0, get_re_lmax),
    Case = case
  )
}

# -----------------------------------------------------------------------------
# Analytical landscape
# -----------------------------------------------------------------------------

grid <- expand_grid(
  T1 = seq(0, 2, length.out = 100),
  T0 = seq(-2, 2, length.out = 100)
) %>%
  mutate(ReLambdaMax = map2_dbl(T1, T0, get_re_lmax))

# -----------------------------------------------------------------------------
# Simulations
# -----------------------------------------------------------------------------

df_all <- bind_rows(
  simulate_scenario("(b) Full dependency", TRUE, TRUE, TRUE),
  simulate_scenario("(c) Feeding-coefficient\ndependency", TRUE, FALSE, FALSE),
  simulate_scenario("(d) Shared-biomass\ndependency", FALSE, TRUE, FALSE),
  simulate_scenario("(e) Shared-efficiency\ndependency", FALSE, FALSE, TRUE),
  simulate_scenario("(f) No dependency", FALSE, FALSE, FALSE)
) %>%
  mutate(Case = factor(Case, levels = case_levels))

x_limits_sim <- range(df_all$T1)
y_limits_sim <- range(df_all$T0)


# -----------------------------------------------------------------------------
# statistical analysis
# -----------------------------------------------------------------------------

# One-way ANOVA
mod_igp <- aov(ReLambdaMax ~ Case, data = df_all)
summary(mod_igp)

# Estimated means and Tukey-adjusted pairwise comparisons
emm_igp <- emmeans(mod_igp, ~ Case)
pairs(emm_igp, adjust = "tukey")

# Letters for the comparison plot
letters_df <- multcomp::cld(
  emm_igp,
  adjust = "tukey",
  Letters = letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group))

# -----------------------------------------------------------------------------
# Panel a: analytical landscape
# -----------------------------------------------------------------------------
base_theme <- theme(
  axis.title = element_text(size = fs, colour = "black"),
  axis.text = element_text(size = fss, colour = "black"),
  plot.title = element_text(size = fss, face = "bold"),
  plot.margin = margin(3, 3, 3, 3, "pt"),
  panel.grid = element_blank()
)

p_analytical <- ggplot(grid, aes(T1, T0, fill = ReLambdaMax)) +
  geom_raster(interpolate = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey35") +
  scale_fill_viridis_c(
    option = "H",
    limits = c(0, 1.4),
    name = expression(Re(lambda[max]))
  ) +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  labs(
    x = expression(T[1]),
    y = expression(T[0]),
    title = "(a) Analytical landscape"
  ) +
  theme_minimal(base_size = fs) +
  base_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = fss),
    legend.text = element_text(size = fss)
  )

# -----------------------------------------------------------------------------
# Ratio panel
# -----------------------------------------------------------------------------

grid_ratio <- grid %>%
  filter(T1 > 0) %>%
  mutate(
    ratio = abs(T0) / abs(T1),
    sign_T0 = if_else(T0 >= 0, "positive", "negative")
  ) %>%
  filter(is.finite(ratio))

p_ratio <- ggplot(grid_ratio, aes(ratio, ReLambdaMax, colour = sign_T0)) +
  geom_point(alpha = 0.2, size = 0.4) +
  scale_colour_manual(
    values = c(positive = "darkblue", negative = "darkred"),
    labels = c(
      positive = expression(T[0] >= 0),
      negative = expression(T[0] < 0)
    ),
    name = expression(sign(T[0]))
  ) +
  scale_x_log10(labels = label_number()) +
  labs(
    x = expression("|" * T[0] * "| / |" * T[1] * "|"),
    y = expression(Re(lambda[max]))
  ) +
  theme_test(base_size = fs) +
  base_theme +
  theme(
    legend.position = c(0.4, 0.96),
    legend.justification = c(1, 1),
    legend.title = element_text(size = fss - 1),
    legend.text = element_text(size = fss - 1),
    legend.background = element_rect(
      fill = alpha("white", 0.75),
      colour = NA
    )
  )

# -----------------------------------------------------------------------------
# Simulated panels
# -----------------------------------------------------------------------------

plot_scenario <- function(case, show_legend = FALSE) {
  df <- filter(df_all, Case == case)
  
  ggplot(df, aes(T1, T0, colour = ReLambdaMax)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = "grey35") +
    geom_point(alpha = 0.45, size = 1) +
    scale_colour_viridis_c(
      option = "H",
      limits = c(0, 0.4),
      breaks = c(0, 0.2, 0.4),
      oob = squish,
      name = expression(Re(lambda[max])),
      guide = if (show_legend) {
        guide_colourbar(
          direction = "horizontal",
          barwidth = unit(3, "cm"),
          barheight = unit(0.25, "cm"),
          title.position = "top"
        )
      } else "none"
    ) +
    coord_cartesian(xlim = x_limits_sim, ylim = y_limits_sim) +
    labs(
      x = expression(T[1]),
      y = expression(T[0]),
      title = case
    ) +
    theme_test(base_size = fs) +
    base_theme +
    theme(
      legend.position = if (show_legend) c(0.7, 0.2) else "none",
      legend.justification = c(0.5, 0.5),
      legend.title = element_text(size = fss - 2),
      legend.text = element_text(size = fss - 3),
      legend.background = element_rect(
        fill = alpha("white", 0.78),
        colour = NA
      )
    )
}

p_full <- plot_scenario("(b) Full dependency")
p_c <- plot_scenario("(c) Feeding-coefficient\ndependency")
p_B <- plot_scenario("(d) Shared-biomass\ndependency")
p_e <- plot_scenario("(e) Shared-efficiency\ndependency")
p_no <- plot_scenario("(f) No dependency", TRUE)

# -----------------------------------------------------------------------------
# Panel g: scenario comparison
# -----------------------------------------------------------------------------

df_box <- df_all %>%
  mutate(Case = factor(Case, levels = case_levels))

p_box <- ggplot(df_box, aes(Case, ReLambdaMax, fill = Case)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  scale_x_discrete(labels = scenario_labs_short) +
  scale_y_continuous(
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, 0.1)
  ) +
  labs(
    x = NULL,
    y = expression(Re(lambda[max])),
    title = "(g) Scenario comparison"
  ) +
  theme_test(base_size = fs) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# -----------------------------------------------------------------------------
# Panel g: scenario comparison
# -----------------------------------------------------------------------------

df_box <- df_all %>%
  mutate(Case = factor(Case, levels = case_levels))

p_box <- ggplot(df_box, aes(Case, ReLambdaMax, fill = Case)) +
  # geom_jitter(
  #   aes(colour = Case),
  #   width = 0.15, height = 0,
  #   alpha = 0.08, size = 0.4
  # ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.7) +
  geom_text(
    data = letters_df,
    aes(x = Case, y = 0.19, label = .group),
    inherit.aes = FALSE,
    size = 5.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  scale_colour_manual(values = scenario_cols, guide = "none") +
  scale_x_discrete(labels = scenario_labs_short) +
  scale_y_continuous(breaks = seq(0, 0.2, 0.1)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(
    x = NULL,
    y = expression(Re(lambda[max])),
    title = "(g) Scenario comparison"
  ) +
  theme_test(base_size = fs) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_box

# -----------------------------------------------------------------------------
# Combine and save
# -----------------------------------------------------------------------------

p_blank <- ggplot() + theme_void()

row1 <- plot_grid(
  p_blank,
  plot_grid(
    p_analytical,
    p_blank,
    p_ratio,
    ncol = 3,
    rel_widths = c(1.1, 0.1, 1.1)
  ),
  p_blank,
  ncol = 3,
  rel_widths = c(0.05, 1, 0.05)
)

p_all <- plot_grid(
  row1,
  plot_grid(p_full, p_c, p_B, ncol = 3),
  plot_grid(p_e, p_no, p_box, ncol = 3),
  ncol = 1,
  rel_heights = c(1.1, 1, 1)
)

p_all

ggsave(here("figures", "3_species_module_dependency.jpeg"),
       p_all, width = 14, height = 13, dpi = 600)
ggsave(here("figures", "3_species_module_dependency.pdf"),
       p_all, width = 14, height = 13)

