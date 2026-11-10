# =============================================================================
# Lake Constance whole-network dependency randomization
# =============================================================================

rm(list = ls())

library(tidyverse)
library(here)
library(colorspace)
library(patchwork)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)

select <- dplyr::select
filter <- dplyr::filter

set.seed(2026)

species <- c("Phyt", "Bac", "HF", "Cil", "Rot", "HerbCru", "CarnCru", "Fish")
n_rep <- 1000

scenario_levels <- c(
  "full_dependency", "c_dependency", "B_dependency",
  "e_dependency", "no_dependency"
)

scenario_labs <- c(
  full_dependency = "Full\ndependency",
  c_dependency    = "Interaction-coefficient\ndependency",
  B_dependency    = "Shared-biomass\ndependency",
  e_dependency    = "Shared-efficiency\ndependency",
  no_dependency   = "No\ndependency"
)

scenario_cols <- c(
  full_dependency = "#009E73",  
  c_dependency    = "#CC79A7",  
  B_dependency    = "#E69F00", 
  e_dependency    = "#D55E00",  
  no_dependency   = "#0072B2"   
)

# -----------------------------------------------------------------------------
# 1. Read biomass and efficiency
# -----------------------------------------------------------------------------

bio_eff <- read_csv(
  here("data", "LC_fws", "biomass_stability_efficiency.csv"),
  show_col_types = FALSE
) %>%
  mutate(web_id = paste(year, phase, sep = "."))

B_list <- bio_eff %>%
  transmute(
    web_id, Phyt, Bac, HF, Cil, Rot, HerbCru, CarnCru,
    Fish = 752
  ) %>%
  split(.$web_id) %>%
  map(~ set_names(as.numeric(.x[1, species]), species))

E_list <- bio_eff %>%
  transmute(
    web_id,
    Phyt = 0,
    Bac = 0,
    HF = HFEfficiency,
    Cil = CilEfficiency,
    Rot = RotEfficiency,
    HerbCru = HerbCruEff,
    CarnCru = CarnCruEff,
    Fish = FishEff
  ) %>%
  split(.$web_id) %>%
  map(~ set_names(as.numeric(.x[1, species]), species))

empirical_metadata <- bio_eff %>%
  transmute(
    web_id, year, phase,
    Re_lambda_file = `Real(lambda)`,
    Lwmax,
    LWmax_3 = `LWmax-3`
  )

# -----------------------------------------------------------------------------
# 2. Read fluxes and calculate feeding coefficients
# -----------------------------------------------------------------------------

flux_df <- read_csv(
  here("data", "LC_fws", "flux_IS.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    web_id = paste(Year, Phase, sep = "."),
    Prey = factor(Prey, levels = 1:8, labels = species),
    Predator = factor(Predator, levels = 1:8, labels = species),
    Flux = replace_na(`Flux mg C m-2 d-1`, 0)
  )

F_list <- split(flux_df, flux_df$web_id) %>%
  map(~ as.matrix(xtabs(
    Flux ~ Prey + Predator,
    data = .x,
    drop.unused.levels = FALSE
  )))

web_ids <- Reduce(intersect, list(names(B_list), names(E_list), names(F_list)))

B_list <- B_list[web_ids]
E_list <- E_list[web_ids]
F_list <- F_list[web_ids]

C_list <- map2(F_list, B_list, \(F, B) {
  C <- F / outer(B, B)
  C[!is.finite(C)] <- 0
  diag(C) <- 0
  dimnames(C) <- list(species, species)
  C
})

# -----------------------------------------------------------------------------
# 3. Jacobian functions
# -----------------------------------------------------------------------------

build_jacobian <- function(C, B, E) {
  J <- matrix(0, length(B), length(B), dimnames = list(names(B), names(B)))
  links <- which(C > 0, arr.ind = TRUE)
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]
    j <- links[k, 2]
    
    J[i, j] <- -C[i, j] * B[i]
    J[j, i] <-  E[j] * C[i, j] * B[j]
  }
  
  J
}

get_Re_lmax <- function(J) max(Re(eigen(J, only.values = TRUE)$values))

# -----------------------------------------------------------------------------
# 4. randomization
# -----------------------------------------------------------------------------

randomize_jacobian <- function(C, B, E, scenario) {
  n <- length(B)
  links <- which(C > 0, arr.ind = TRUE)
  
  c_pool <- as.numeric(C[C > 0])
  B_pool <- as.numeric(B)
  e_pool <- as.numeric(E[E > 0])
  
  J <- matrix(0, n, n, dimnames = list(names(B), names(B)))
  
  share_c <- scenario %in% c("full_dependency", "c_dependency")
  share_B <- scenario %in% c("full_dependency", "B_dependency")
  share_e <- scenario %in% c("full_dependency", "e_dependency")
  
  if (share_B) B_new <- sample(B_pool, n, replace = TRUE)
  
  if (share_e) {
    E_new <- numeric(n)
    consumers <- which(E > 0)
    E_new[consumers] <- sample(e_pool, length(consumers), replace = TRUE)
  }
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]
    j <- links[k, 2]
    
    if (share_c) {
      c_ij <- sample(c_pool, 1, replace = TRUE)
      c_top <- c_ij
      c_bottom <- c_ij
    } else {
      c_top <- sample(c_pool, 1, replace = TRUE)
      c_bottom <- sample(c_pool, 1, replace = TRUE)
    }
    
    if (share_B) {
      B_top <- B_new[i]
      B_bottom <- B_new[j]
    } else {
      B_top <- sample(B_pool, 1, replace = TRUE)
      B_bottom <- sample(B_pool, 1, replace = TRUE)
    }
    
    if (share_e) {
      e_bottom <- E_new[j]
    } else {
      e_bottom <- sample(e_pool, 1, replace = TRUE)
    }
    
    J[i, j] <- J[i, j] - c_top * B_top
    J[j, i] <- J[j, i] + e_bottom * c_bottom * B_bottom
  }
  
  J
}

# -----------------------------------------------------------------------------
# 5. Empirical and randomized stability
# -----------------------------------------------------------------------------

empirical_results <- map_dfr(web_ids, \(id) {
  tibble(
    web_id = id,
    Re_emp = get_Re_lmax(build_jacobian(
      C_list[[id]], B_list[[id]], E_list[[id]]
    ))
  )
}) %>%
  left_join(empirical_metadata, by = "web_id")

randomization_results <- map_dfr(web_ids, \(id) {
  message("Running ", id)
  
  crossing(
    scenario = scenario_levels,
    rep = seq_len(n_rep)
  ) %>%
    mutate(
      web_id = id,
      Re_lmax = map_dbl(
        scenario,
        \(s) get_Re_lmax(randomize_jacobian(
          C_list[[id]], B_list[[id]], E_list[[id]], s
        ))
      )
    )
})

# -----------------------------------------------------------------------------
# 6. Scenario summaries and effect sizes
# -----------------------------------------------------------------------------
scenario_summary <- randomization_results %>%
  group_by(web_id, scenario) %>%
  summarise(
    mean_Re = mean(Re_lmax),
    sd_Re = sd(Re_lmax),
    .groups = "drop"
  ) %>%
  left_join(select(empirical_results, web_id, Re_emp), by = "web_id") %>%
  mutate(
    raw_effect = mean_Re - Re_emp,
    standardized_effect = raw_effect / sd_Re,
    scenario = factor(scenario, levels = scenario_levels)
  )

# Positive effect:
# randomized Re(lambda_max) > empirical Re(lambda_max)

check_df <- scenario_summary %>%
  mutate(
    difference = mean_Re - Re_emp,
    abs_difference = abs(difference),
    abs_SES = abs(standardized_effect)
  )

# How many food webs have randomised means above the empirical value?
count_summary <- check_df %>%
  group_by(scenario) %>%
  summarise(
    n_webs = n(),
    n_above_empirical = sum(difference > 0),
    n_below_empirical = sum(difference < 0),
    n_equal_empirical = sum(difference == 0),
    .groups = "drop"
  )

print(count_summary)

# -----------------------------------------------------------------------------
# 7. Statistical analysis
# -----------------------------------------------------------------------------

mod_lc <- lmerTest::lmer(
  standardized_effect ~ scenario + (1 | web_id),
  data = scenario_summary
)

# ANOVA-style test of scenario differences
anova_result <- anova(mod_lc, ddf = "Satterthwaite")
print(anova_result)

# Estimated means and Tukey-adjusted pairwise comparisons
emm_lc <- emmeans(
  mod_lc, ~ scenario,
  lmer.df = "satterthwaite"
)

pairwise_tests <- pairs(emm_lc, adjust = "tukey")
print(pairwise_tests)

# Letters for the SES plot
letters_df <- multcomp::cld(
  emm_lc,
  adjust = "tukey",
  alpha = 0.05,
  Letters = letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group)) %>%
  select(scenario, .group)


# -----------------------------------------------------------------------------
# 8. Plotting
# -----------------------------------------------------------------------------

time_index <- empirical_results %>%
  distinct(web_id, year, phase, Re_emp) %>%
  group_by(year) %>%
  arrange(phase, .by_group = TRUE) %>%
  mutate(
    phase_order = row_number(),
    year_time = year + (phase_order - (n() + 1) / 2) / (n() + 1)
  ) %>%
  ungroup()

time_df <- scenario_summary %>%
  left_join(time_index, by = "web_id")

fs <- 18

p_time <- ggplot(
  time_df,
  aes(year_time, mean_Re, colour = scenario, group = scenario)
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1) +
  geom_line(
    data = time_index,
    aes(year_time, Re_emp, group = 1),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.9
  ) +
  geom_point(
    data = time_index,
    aes(year_time, Re_emp),
    inherit.aes = FALSE,
    colour = "black",
    size = 1
  ) +
  scale_colour_manual(values = scenario_cols, labels = scenario_labs) +
  scale_x_continuous(
    breaks = sort(unique(time_index$year)),
    labels = sort(unique(time_index$year))
  ) +
  labs(
    x = "Year",
    y = expression(Re(lambda[max])),
    colour = NULL
  ) +
  theme_classic(base_size = fs) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = fs - 2),
    axis.text = element_text(colour = "black")
  )

# Letter positions
y_gap <- 0.09 * max(
  1,
  diff(range(scenario_summary$standardized_effect, na.rm = TRUE))
)

letter_positions <- scenario_summary %>%
  group_by(scenario) %>%
  summarise(
    label_y = max(standardized_effect, na.rm = TRUE) + y_gap,
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = "scenario")

p_ses <- ggplot(
  scenario_summary,
  aes(scenario, standardized_effect, fill = scenario, colour = scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.68
  ) +
  geom_jitter(
    width = 0.12, height = 0,
    size = 1.4, alpha = 0.7
  ) +
  geom_text(
    data = letter_positions,
    aes(x = scenario, y = label_y, label = .group),
    inherit.aes = FALSE,
    colour = "black",
    size = 5.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = scenario_cols) +
  scale_colour_manual(values = scenario_cols) +
  scale_x_discrete(labels = scenario_labs) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  guides(fill = "none", colour = "none") +
  labs(x = NULL, y = "SES") +
  theme_classic(base_size = fs) +
  theme(
    panel.border = element_rect(
      colour = "black", fill = NA, linewidth = 0.45
    ),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(colour = "black")
  )

final_plot <- p_time / p_ses +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
final_plot

ggsave(here("figures", "Figure3.jpeg"), 
       final_plot, width = 12, height = 10, dpi = 600)
ggsave(here("figures", "Figure3.pdf"), 
       final_plot, width = 12, height = 10)

#write_csv(scenario_summary, here("data", "LC_fws", "dependency_randomization_summary.csv"))
