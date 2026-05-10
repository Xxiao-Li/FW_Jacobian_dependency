#============================================================
# 2. Whole-network dependency randomisations for Lake Constance
#============================================================

rm(list = ls())

library(tidyverse)
library(here)
library(purrr)
library(scico)
library(lme4)
library(lmerTest)
library(emmeans)
library(scales)
library(colorspace)
library(cowplot)
library(patchwork)
library(magick)
library(cowplot)

set.seed(2026)

#------------------------------------------------------------
# 1. Basic settings
#------------------------------------------------------------

# Species order:
# 1 = Phyt, 2 = Bac, 3 = HF, 4 = Cil, 5 = Rot,
# 6 = HerbCru, 7 = CarnCru, 8 = Fish

species <- c(
  "Phyt", "Bac", "HF", "Cil",
  "Rot", "HerbCru", "CarnCru", "Fish"
)

n_rep <- 1000

#------------------------------------------------------------
# 2. Read biomass, efficiency, empirical stability
#------------------------------------------------------------

bio_eff_df <- read_csv(here("data", "biomass_stability_efficiency.csv"))

# Biomass
spe_biomass <- bio_eff_df %>%
  select(year, phase, Phyt, Bac, HF, Cil, Rot, HerbCru, CarnCru) %>%
  mutate(Fish = 752.00)

B_list <- spe_biomass %>%
  mutate(web_id = paste(year, phase, sep = ".")) %>%
  split(.$web_id) %>%
  lapply(function(df) {
    B <- as.numeric(df[1, species])
    names(B) <- species
    B
  })

# Efficiency
# Basal species Phyt and Bac are assigned efficiency = 0
spe_efficiency <- bio_eff_df %>%
  select(
    year, phase,
    HFEfficiency, CilEfficiency, RotEfficiency,
    HerbCruEff, CarnCruEff, FishEff
  )

E_list <- spe_efficiency %>%
  mutate(web_id = paste(year, phase, sep = ".")) %>%
  split(.$web_id) %>%
  lapply(function(df) {
    E <- c(
      Phyt    = 0,
      Bac     = 0,
      HF      = df$HFEfficiency[1],
      Cil     = df$CilEfficiency[1],
      Rot     = df$RotEfficiency[1],
      HerbCru = df$HerbCruEff[1],
      CarnCru = df$CarnCruEff[1],
      Fish    = df$FishEff[1]
    )
    E
  })

# Empirical stability and loop-weight values if needed for comparison
fw_stability_LWmax <- bio_eff_df %>%
  select(year, phase, "Real(lambda)", Lwmax, "LWmax-3") %>%
  rename(
    Re_lambda_emp = "Real(lambda)",
    LWmax_3       = "LWmax-3"
  ) %>%
  mutate(web_id = paste(year, phase, sep = "."))

#------------------------------------------------------------
# 3. Read flux data and convert to flux matrices
#------------------------------------------------------------

flux_df <- read_csv(here("data", "flux_IS.csv")) %>%
  select(Year, Phase, Prey, Predator, "Flux mg C m-2 d-1") %>%
  rename(Flux = "Flux mg C m-2 d-1") %>%
  mutate(
    Flux = if_else(is.na(Flux), 0, Flux),
    Prey = factor(Prey, levels = 1:8, labels = species),
    Predator = factor(Predator, levels = 1:8, labels = species),
    web_id = paste(Year, Phase, sep = ".")
  )

flux_by_web <- split(flux_df, flux_df$web_id)

flux_mat_list <- lapply(flux_by_web, function(df) {
  m <- xtabs(Flux ~ Prey + Predator, data = df, drop = FALSE)
  m <- as.matrix(m)
  dimnames(m) <- list(
    Prey = species,
    Predator = species
  )
  m
})

# Keep only webs present in all datasets
web_ids <- Reduce(
  intersect,
  list(names(B_list), names(E_list), names(flux_mat_list))
)

B_list <- B_list[web_ids]
E_list <- E_list[web_ids]
flux_mat_list <- flux_mat_list[web_ids]

saveRDS(B_list, here("data", "B_list.rds"))
saveRDS(E_list, here("data", "E_list.rds"))
saveRDS(flux_mat_list, here("data", "flux_mat_list.rds"))

#------------------------------------------------------------
# 4. Calculate feeding coefficients c_ij
#------------------------------------------------------------
# c_ij = F_ij / (B_i * B_j)
# Rows are resources/prey i, columns are consumers/predators j

cij_list <- imap(flux_mat_list, function(F, web_id) {
  B <- B_list[[web_id]]
  
  den <- B %o% B
  C <- F / den
  
  C[!is.finite(C)] <- 0
  diag(C) <- 0
  
  dimnames(C) <- dimnames(F)
  
  C
})

#------------------------------------------------------------
# 5. Build Jacobian matrix
#------------------------------------------------------------
# For a trophic link where consumer j feeds on resource i:
# top-down effect:  J_ij = -c_ij * B_i
# bottom-up effect: J_ji =  e_j * c_ij * B_j
#
# This function uses all non-zero entries in C, not only the upper triangle.

build_jacobian <- function(C, B, E) {
  n <- length(B)
  
  stopifnot(
    is.matrix(C),
    all(dim(C) == c(n, n)),
    length(E) == n
  )
  
  J <- matrix(0, n, n)
  rownames(J) <- names(B)
  colnames(J) <- names(B)
  
  links <- which(C > 0, arr.ind = TRUE)
  
  if (nrow(links) > 0) {
    for (k in seq_len(nrow(links))) {
      i <- links[k, 1]  # prey/resource row
      j <- links[k, 2]  # predator/consumer column
      
      J[i, j] <- J[i, j] - C[i, j] * B[i]
      J[j, i] <- J[j, i] + E[j] * C[i, j] * B[j]
    }
  }
  
  diag(J) <- 0
  
  J
}

get_Re_lmax <- function(J) {
  ev <- eigen(J, only.values = TRUE)$values
  max(Re(ev))
}

compute_Re_lmax_from_params <- function(C, B, E) {
  J <- build_jacobian(C, B, E)
  get_Re_lmax(J)
}

#------------------------------------------------------------
# 6. Empirical Jacobian stability
#------------------------------------------------------------

empirical_results <- map_dfr(web_ids, function(web_id) {
  C <- cij_list[[web_id]]
  B <- B_list[[web_id]]
  E <- E_list[[web_id]]
  
  tibble(
    web_id = web_id,
    scenario = "Empirical",
    rep = NA_integer_,
    Re_lmax = compute_Re_lmax_from_params(C, B, E)
  )
}) %>%
  left_join(
    fw_stability_LWmax %>%
      select(web_id, year, phase, Re_lambda_emp, Lwmax, LWmax_3),
    by = "web_id"
  )

# check: compare computed empirical Re_lmax with existing values
empirical_check <- empirical_results %>%
  mutate(diff_from_file = Re_lmax - Re_lambda_emp)

print(summary(empirical_check$diff_from_file))

#------------------------------------------------------------
# 7. Helper functions for randomisation
#------------------------------------------------------------

# Extract empirical parameter pools for one web
get_parameter_pools <- function(C, B, E) {
  c_pool <- as.numeric(C[C > 0])
  B_pool <- as.numeric(B)
  
  # Here we use E > 0 because Phyt and Bac have efficiency = 0
  e_pool <- as.numeric(E[E > 0])
  
  list(
    c_pool = c_pool,
    B_pool = B_pool,
    e_pool = e_pool
  )
}

get_links <- function(C) {
  links <- which(C > 0, arr.ind = TRUE)
  
  tibble(
    prey = links[, 1],
    predator = links[, 2]
  )
}

empty_params <- function(n, species) {
  C_new <- matrix(0, n, n)
  dimnames(C_new) <- list(species, species)
  
  B_new <- rep(NA_real_, n)
  names(B_new) <- species
  
  E_new <- rep(NA_real_, n)
  names(E_new) <- species
  
  list(C = C_new, B = B_new, E = E_new)
}

#------------------------------------------------------------
# 8. Randomisation scenarios
#------------------------------------------------------------
# Scenario names:
# full_dependency
# c_dependency
# B_dependency
# e_dependency
# no_dependency
#
# In all scenarios, empirical topology and signs are preserved.
# Values are drawn from the empirical pools of each food web.

randomise_params <- function(C, B, E, scenario) {
  n <- length(B)
  sp <- names(B)
  
  pools <- get_parameter_pools(C, B, E)
  links <- get_links(C)
  
  c_pool <- pools$c_pool
  B_pool <- pools$B_pool
  e_pool <- pools$e_pool
  
  if (length(c_pool) == 0) stop("c_pool is empty.")
  if (length(B_pool) == 0) stop("B_pool is empty.")
  if (length(e_pool) == 0) stop("e_pool is empty.")
  
  #----------------------------------------------------------
  # Full structural dependency
  #----------------------------------------------------------
  # One B per species, one E per consumer species,
  # one C per trophic link. Reused consistently in J.
  if (scenario == "full_dependency") {
    C_new <- matrix(0, n, n, dimnames = list(sp, sp))
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      C_new[i, j] <- sample(c_pool, 1, replace = TRUE)
    }
    
    B_new <- sample(B_pool, n, replace = TRUE)
    names(B_new) <- sp
    
    E_new <- rep(0, n)
    names(E_new) <- sp
    
    consumer_idx <- which(E > 0)
    E_new[consumer_idx] <- sample(e_pool, length(consumer_idx), replace = TRUE)
    
    return(list(C = C_new, B = B_new, E = E_new))
  }
  
  #----------------------------------------------------------
  # Consumer-resource interaction dependency only
  #----------------------------------------------------------
  # C is reused within each predator-prey pair.
  # B and E are not represented by one shared vector.
  # Instead, we directly build the Jacobian.
  if (scenario == "c_dependency") {
    J <- matrix(0, n, n, dimnames = list(sp, sp))
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      c_ij <- sample(c_pool, 1, replace = TRUE)
      
      B_topdown <- sample(B_pool, 1, replace = TRUE)
      B_bottomup <- sample(B_pool, 1, replace = TRUE)
      e_bottomup <- sample(e_pool, 1, replace = TRUE)
      
      J[i, j] <- J[i, j] - c_ij * B_topdown
      J[j, i] <- J[j, i] + e_bottomup * c_ij * B_bottomup
    }
    
    diag(J) <- 0
    return(list(J = J))
  }
  
  #----------------------------------------------------------
  # Shared biomass dependency only
  #----------------------------------------------------------
  # One B per species is reused.
  # Each occurrence of C and E is independently resampled.
  if (scenario == "B_dependency") {
    B_new <- sample(B_pool, n, replace = TRUE)
    names(B_new) <- sp
    
    J <- matrix(0, n, n, dimnames = list(sp, sp))
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      c_topdown <- sample(c_pool, 1, replace = TRUE)
      c_bottomup <- sample(c_pool, 1, replace = TRUE)
      e_bottomup <- sample(e_pool, 1, replace = TRUE)
      
      J[i, j] <- J[i, j] - c_topdown * B_new[i]
      J[j, i] <- J[j, i] + e_bottomup * c_bottomup * B_new[j]
    }
    
    diag(J) <- 0
    return(list(J = J))
  }
  
  #----------------------------------------------------------
  # Shared consumer efficiency dependency only
  #----------------------------------------------------------
  # One E per consumer is reused.
  # Each occurrence of C and B is independently resampled.
  if (scenario == "e_dependency") {
    E_new <- rep(0, n)
    names(E_new) <- sp
    
    consumer_idx <- which(E > 0)
    E_new[consumer_idx] <- sample(e_pool, length(consumer_idx), replace = TRUE)
    
    J <- matrix(0, n, n, dimnames = list(sp, sp))
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      c_topdown <- sample(c_pool, 1, replace = TRUE)
      B_topdown <- sample(B_pool, 1, replace = TRUE)
      
      c_bottomup <- sample(c_pool, 1, replace = TRUE)
      B_bottomup <- sample(B_pool, 1, replace = TRUE)
      
      J[i, j] <- J[i, j] - c_topdown * B_topdown
      J[j, i] <- J[j, i] + E_new[j] * c_bottomup * B_bottomup
    }
    
    diag(J) <- 0
    return(list(J = J))
  }
  
  #----------------------------------------------------------
  # No dependency
  #----------------------------------------------------------
  # Every occurrence of C, B, and E is independently sampled.
  if (scenario == "no_dependency") {
    J <- matrix(0, n, n, dimnames = list(sp, sp))
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      c_topdown <- sample(c_pool, 1, replace = TRUE)
      B_topdown <- sample(B_pool, 1, replace = TRUE)
      
      c_bottomup <- sample(c_pool, 1, replace = TRUE)
      B_bottomup <- sample(B_pool, 1, replace = TRUE)
      e_bottomup <- sample(e_pool, 1, replace = TRUE)
      
      J[i, j] <- J[i, j] - c_topdown * B_topdown
      J[j, i] <- J[j, i] + e_bottomup * c_bottomup * B_bottomup
    }
    
    diag(J) <- 0
    return(list(J = J))
  }
  
  stop("Unknown scenario: ", scenario)
}

compute_randomised_Re_lmax <- function(C, B, E, scenario) {
  rand <- randomise_params(C, B, E, scenario)
  
  if (!is.null(rand$J)) {
    return(get_Re_lmax(rand$J))
  }
  
  J <- build_jacobian(rand$C, rand$B, rand$E)
  get_Re_lmax(J)
}

#------------------------------------------------------------
# 9. Run randomisations
#------------------------------------------------------------

scenarios <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

randomisation_results <- map_dfr(web_ids, function(web_id) {
  message("Running web: ", web_id)
  
  C <- cij_list[[web_id]]
  B <- B_list[[web_id]]
  E <- E_list[[web_id]]
  
  map_dfr(scenarios, function(scen) {
    message("  Scenario: ", scen)
    
    tibble(rep = seq_len(n_rep)) %>%
      mutate(
        web_id = web_id,
        scenario = scen,
        Re_lmax = map_dbl(
          rep,
          ~ compute_randomised_Re_lmax(C, B, E, scen)
        )
      )
  })
})


#------------------------------------------------------------
# 10. Compare randomisations with empirical food webs
#------------------------------------------------------------

# Mean and SD of randomised stability values per web and scenario
scenario_summary_by_web <- randomisation_results %>%
  group_by(web_id, scenario) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    sd_Re   = sd(Re_lmax, na.rm = TRUE),
    median_Re = median(Re_lmax, na.rm = TRUE),
    q025_Re = quantile(Re_lmax, 0.025, na.rm = TRUE),
    q975_Re = quantile(Re_lmax, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Add empirical stability
effect_size_results <- scenario_summary_by_web %>%
  left_join(
    empirical_results %>%
      select(web_id, Re_emp = Re_lmax),
    by = "web_id"
  ) %>%
  mutate(
    raw_effect = mean_Re - Re_emp,
    standardised_effect = (Re_emp - mean_Re) / sd_Re
  )

# Interpretation:
# raw_effect > 0:
#   randomised matrices have higher Re_lmax than empirical food web.
#
# standardised_effect < 0:
#   empirical food web has lower Re_lmax than null expectation.

# Test whether raw effects differ from zero across food webs
effect_tests <- effect_size_results %>%
  group_by(scenario) %>%
  summarise(
    n_webs = n(),
    mean_raw_effect = mean(raw_effect, na.rm = TRUE),
    median_raw_effect = median(raw_effect, na.rm = TRUE),
    sd_raw_effect = sd(raw_effect, na.rm = TRUE),
    t_test_p = t.test(raw_effect, mu = 0)$p.value,
    wilcox_p = wilcox.test(raw_effect, mu = 0)$p.value,
    mean_standardised_effect = mean(standardised_effect, na.rm = TRUE),
    median_standardised_effect = median(standardised_effect, na.rm = TRUE),
    .groups = "drop"
  )


#------------------------------------------------------------
# 11. Compare randomisation scenarios with each other
#------------------------------------------------------------

# Use one value per food web per scenario to avoid treating
# 1000 randomisations as independent biological replicates.

scenario_means <- randomisation_results %>%
  group_by(web_id, scenario) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "no_dependency",
        "c_dependency",
        "B_dependency",
        "e_dependency",
        "full_dependency"
      )
    )
  )

# Repeated-measures ANOVA
anova_fit <- aov(
  mean_Re ~ scenario + Error(web_id / scenario),
  data = scenario_means
)

anova_summary <- summary(anova_fit)
print(anova_summary)

# Linear mixed-effects model
lmm_fit <- lmer(
  mean_Re ~ scenario + (1 | web_id),
  data = scenario_means
)

print(summary(lmm_fit))
print(anova(lmm_fit))

# Post-hoc pairwise comparisons among scenarios
emm <- emmeans(lmm_fit, ~ scenario)
pairwise_scenario_comparisons <- pairs(emm, adjust = "tukey")

print(pairwise_scenario_comparisons)

pairwise_df <- as.data.frame(pairwise_scenario_comparisons)


# Specific contrasts against no dependency
contrasts_vs_no <- contrast(
  emm,
  method = list(
    "c_dependency - no_dependency" =
      c(-1, 1, 0, 0, 0),
    "B_dependency - no_dependency" =
      c(-1, 0, 1, 0, 0),
    "e_dependency - no_dependency" =
      c(-1, 0, 0, 1, 0),
    "full_dependency - no_dependency" =
      c(-1, 0, 0, 0, 1)
  ),
  adjust = "holm"
)

contrasts_vs_no_df <- as.data.frame(contrasts_vs_no)

#------------------------------------------------------------
# 12. Dependency-effect contrasts
#------------------------------------------------------------
# Positive values mean the dependency scenario has lower Re_lmax
# than the no-dependency scenario.

dependency_contrasts <- scenario_means %>%
  pivot_wider(
    names_from = scenario,
    values_from = mean_Re
  ) %>%
  mutate(
    effect_c_dependency =
      no_dependency - c_dependency,
    effect_B_dependency =
      no_dependency - B_dependency,
    effect_e_dependency =
      no_dependency - e_dependency,
    effect_full_dependency =
      no_dependency - full_dependency
  )


dependency_contrast_tests <- dependency_contrasts %>%
  summarise(
    c_dep_mean = mean(effect_c_dependency, na.rm = TRUE),
    c_dep_t_p = t.test(effect_c_dependency, mu = 0)$p.value,
    c_dep_wilcox_p = wilcox.test(effect_c_dependency, mu = 0)$p.value,
    
    B_dep_mean = mean(effect_B_dependency, na.rm = TRUE),
    B_dep_t_p = t.test(effect_B_dependency, mu = 0)$p.value,
    B_dep_wilcox_p = wilcox.test(effect_B_dependency, mu = 0)$p.value,
    
    e_dep_mean = mean(effect_e_dependency, na.rm = TRUE),
    e_dep_t_p = t.test(effect_e_dependency, mu = 0)$p.value,
    e_dep_wilcox_p = wilcox.test(effect_e_dependency, mu = 0)$p.value,
    
    full_dep_mean = mean(effect_full_dependency, na.rm = TRUE),
    full_dep_t_p = t.test(effect_full_dependency, mu = 0)$p.value,
    full_dep_wilcox_p = wilcox.test(effect_full_dependency, mu = 0)$p.value
  )


#------------------------------------------------------------
# 13. Plots
#------------------------------------------------------------
scenario_levels <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

scenario_labs <- c(
  full_dependency = "Full structural\ndependency",
  c_dependency    = "Consumer-resource\nc dependency",
  B_dependency    = "Shared biomass\nB dependency",
  e_dependency    = "Shared efficiency\ne dependency",
  no_dependency   = "No dependency"
)

scenario_cols <- scico(
  n = 5,
  palette = "batlow",
  begin = 0.1,
  end = 0.9
)

names(scenario_cols) <- scenario_levels

scenario_fills <- lighten(scenario_cols, amount = 0.65)
scenario_fills <- alpha(scenario_fills, 0.18)

# ------------------------------------------------------------------
# Summarise randomisation results
# ------------------------------------------------------------------

scenario_summary_by_web <- randomisation_results %>%
  group_by(web_id, scenario) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    sd_Re   = sd(Re_lmax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  )

# ------------------------------------------------------------------
# Panel a: temporal comparison
# ------------------------------------------------------------------

# First create a unique time index for each seasonal food web
time_index <- empirical_results %>%
  select(web_id, year, phase, Re_emp = Re_lmax) %>%
  distinct() %>%
  group_by(year) %>%
  arrange(year, phase, .by_group = TRUE) %>%
  mutate(
    phase_order = dense_rank(phase),
    n_phase = n_distinct(phase),
    year_time = year + (phase_order - (n_phase + 1) / 2) / (n_phase + 1)
  ) %>%
  ungroup()

# Join with scenario summaries
time_df <- time_index %>%
  left_join(
    scenario_summary_by_web,
    by = "web_id"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  )

# Empirical line data
emp_time_df <- time_index

fs <- 16

p_time <- time_df %>%
  ggplot(aes(x = year_time, y = mean_Re, color = scenario, fill = scenario, group = scenario)) +
  # geom_ribbon(
  #   aes(
  #     ymin = mean_Re - sd_Re,
  #     ymax = mean_Re + sd_Re
  #   ),
  #   color = NA,
  #   alpha = 0.15
  # ) +
  geom_line(linewidth = 0.75, alpha = 0.95) +
  geom_point(size = 1.0, alpha = 0.85) +
  geom_line(data = emp_time_df,
    aes(x = year_time, y = Re_emp, group = 1),
    inherit.aes = FALSE, color = "black", linewidth = 0.9 ) +
  geom_point(data = emp_time_df,
    aes(x = year_time, y = Re_emp),
    inherit.aes = FALSE, color = "black", size = 1.0) +
  scale_x_continuous(
    breaks = sort(unique(time_index$year)),
    labels = sort(unique(time_index$year))
  ) +
  scale_color_manual(
    values = scenario_cols,
    labels = scenario_labs
  ) +
  scale_fill_manual(
    values = scenario_fills,
    guide = "none"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x = "Year",
    y = expression("Re("*lambda[max]*")"),
    color = "Scenario"
  ) +
  theme_bw(base_size = fs) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = fs-2, color = "black"),
    legend.text = element_text(size = fs - 2, color = "black"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_text(size = fs, color = "black"),
    axis.text = element_text(size = fs - 2, color = "black")
  )

p_time
# ------------------------------------------------------------------
# Panel b: standardised effect size
# ------------------------------------------------------------------
p_ses <- effect_size_results %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  ) %>%
  ggplot(aes(x = scenario, y = standardised_effect, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = c(-2, 1), linetype = "dotted", alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.65) +
  geom_point(position = position_jitter(width = 0.1),
    size = 1.5, alpha = 0.35, color = "gray60", fill = "gray85" ) +
  scale_x_discrete(labels = scenario_labs) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  labs(
    x = NULL,
    y = "Standardised effect size"
  ) +
  theme_bw(base_size = fs) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text = element_text(size= fs, color = "black"),
    axis.title = element_text(size = fs, color = "black"),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

# ------------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------------
final_plot <- p_time / p_ses +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = fs)
  )

final_plot

# ------------------------------------------------------------------
# Save figure
# ------------------------------------------------------------------

ggsave(filename = here("figures", "LakeConstance_dependency.png"),
  plot = final_plot, width = 10, height = 10, dpi = 600)

ggsave(filename = here("figures", "LakeConstance_dependency.pdf"),
  plot = final_plot, width = 10, height = 10)


# ------------------------------------------------------------------
### statistics in the text
# ------------------------------------------------------------------
scenario_levels <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

# Empirical Re(lambda_max) summary
emp_stats <- empirical_results %>%
  summarise(
    n_web = n(),
    median_emp = median(Re_lmax, na.rm = TRUE),
    q25_emp = quantile(Re_lmax, 0.25, na.rm = TRUE),
    q75_emp = quantile(Re_lmax, 0.75, na.rm = TRUE),
    mean_emp = mean(Re_lmax, na.rm = TRUE),
    sd_emp = sd(Re_lmax, na.rm = TRUE))

emp_stats

# Randomised Re(lambda_max) summary by scenario
rand_stats <- scenario_summary_by_web %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  ) %>%
  group_by(scenario) %>%
  summarise(
    n_web = n(),
    median_Re = median(mean_Re, na.rm = TRUE),
    q25_Re = quantile(mean_Re, 0.25, na.rm = TRUE),
    q75_Re = quantile(mean_Re, 0.75, na.rm = TRUE),
    mean_Re = mean(mean_Re, na.rm = TRUE),
    sd_Re = sd(mean_Re, na.rm = TRUE),
    .groups = "drop"
  )

rand_stats

# Standardised effect size summary
ses_stats <- effect_size_results %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  ) %>%
  group_by(scenario) %>%
  summarise(
    n_web = n(),
    median_SES = median(standardised_effect, na.rm = TRUE),
    q25_SES = quantile(standardised_effect, 0.25, na.rm = TRUE),
    q75_SES = quantile(standardised_effect, 0.75, na.rm = TRUE),
    mean_SES = mean(standardised_effect, na.rm = TRUE),
    sd_SES = sd(standardised_effect, na.rm = TRUE),
    n_negative = sum(standardised_effect < 0, na.rm = TRUE),
    prop_negative = mean(standardised_effect < 0, na.rm = TRUE),
    n_below_minus2 = sum(standardised_effect < -2, na.rm = TRUE),
    prop_below_minus2 = mean(standardised_effect < -2, na.rm = TRUE),
    .groups = "drop"
  )

ses_stats

# Raw empirical-randomised difference by scenario
# positive values mean randomised Re(lambda_max) > empirical Re(lambda_max)
raw_diff_stats <- effect_size_results %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    raw_diff = mean_Re - Re_emp
  ) %>%
  group_by(scenario) %>%
  summarise(
    n_web = n(),
    median_raw_diff = median(raw_diff, na.rm = TRUE),
    q25_raw_diff = quantile(raw_diff, 0.25, na.rm = TRUE),
    q75_raw_diff = quantile(raw_diff, 0.75, na.rm = TRUE),
    mean_raw_diff = mean(raw_diff, na.rm = TRUE),
    sd_raw_diff = sd(raw_diff, na.rm = TRUE),
    n_emp_lower = sum(raw_diff > 0, na.rm = TRUE),
    prop_emp_lower = mean(raw_diff > 0, na.rm = TRUE),
    .groups = "drop"
  )

raw_diff_stats

# Paired tests comparing every two dependency scenarios
scenario_levels <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

# Wide format: one row per food web, one column per scenario
scenario_wide <- scenario_summary_by_web %>%
  select(web_id, scenario, mean_Re) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = mean_Re
  )

# Overall repeated-measures test across all scenarios
friedman_result <- friedman.test(
  as.matrix(scenario_wide[, scenario_levels])
)

friedman_result

# Pairwise paired Wilcoxon tests among all scenario pairs
scenario_pairs <- combn(scenario_levels, 2, simplify = FALSE)

pairwise_scenario_tests_tidy <- map_dfr(scenario_pairs, function(pair) {
  s1 <- pair[1]
  s2 <- pair[2]
  
  x <- scenario_wide[[s1]]
  y <- scenario_wide[[s2]]
  
  tibble(
    scenario_1 = s1,
    scenario_2 = s2,
    median_difference = median(x - y, na.rm = TRUE),
    n_scenario_1_lower = sum(x < y, na.rm = TRUE),
    n_web = sum(!is.na(x) & !is.na(y)),
    p_value = wilcox.test(x, y, paired = TRUE)$p.value
  )
}) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "holm")
  )

pairwise_scenario_tests_tidy

