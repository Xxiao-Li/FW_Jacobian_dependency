#============================================================
# Effects of dependency in model food webs
#============================================================

rm(list = ls())

library(tidyverse)
library(here)
library(arrow)
library(patchwork)
library(broom)
library(lme4)
library(emmeans)
library(multcomp)
library(scales)
library(lmerTest)
library(multcompView)

select <- dplyr::select
filter <- dplyr::filter

set.seed(2026)

source(here("code/src/internal.R"))

#------------------------------------------------------------
# 1. Settings
#------------------------------------------------------------
out_dir <- here("data/model_fws")
feather_file <- file.path(out_dir, "model_food_webs_dependency_links.feather") # richness 10-40

# TEST baseline stability:
# "Re_emp"     = empirical-style Jacobian, diagonal = 0
# "Re_intdiag" = off-diagonal trophic terms + interaction-derived diagonal only
# "Re_ode"     = full ODE Jacobian, including interaction-derived and background diagonal

# reference_type <- "Re_emp" # To do randomization and to quantify stabilizing effect of full dependency
# reference_type <- "Re_intdiag"
  reference_type <- "Re_ode"

reference_label <- recode(
  reference_type,
  Re_emp     = "Re(emp)",
  Re_intdiag = "Re(int. diag)",
  Re_ode     = "Re(ode)"
)

reference_suffix <- paste0("_", reference_type)

randomisation_rds <- file.path(out_dir, 
                               paste0("model_fw_dependency_randomization", reference_suffix, ".rds"))

fig_jpeg <- here("figures", paste0("model_fws_dependency", reference_suffix, ".jpeg"))
fig_pdf  <- here("figures", paste0("model_fws_dependency", reference_suffix, ".pdf"))

S_label_levels <- c("6_10", "11_15", "16_20", "21_25", "26_30", "31_35", "36_40")

S_label_labs <- c(
  "6_10"  = "6–10",
  "11_15" = "11–15",
  "16_20" = "16–20",
  "21_25" = "21–25",
  "26_30" = "26–30",
  "31_35" = "31–35",
  "36_40" = "36–40"
)

scenario_levels <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

scenario_labs <- c(
  full_dependency = "Full dependency",
  c_dependency    = "Feeding-coefficient\ndependency",
  B_dependency    = "Shared-biomass\ndependency",
  e_dependency    = "Shared-efficiency\ndependency",
  no_dependency   = "No dependency"
)

scenario_cols <- c(
  full_dependency = "#009E73",
  c_dependency    = "#CC79A7",
  B_dependency    = "#E69F00",
  e_dependency    = "#D55E00",
  no_dependency   = "#0072B2"
)

scenario_cols_fill <- scales::alpha(scenario_cols, 0.25)

n_rep <- 1000 # for each scenario and each model food web

#------------------------------------------------------------
# 2. Read Julia output and rebuild webs
#------------------------------------------------------------
modelfw_df <- read_feather(feather_file) %>%
  mutate(
    valid_id = as.character(valid_id),
    S_label = str_remove(as.character(S_label), "^S"),
    S_label = factor(S_label, levels = S_label_levels)
  )

web_objects <- modelfw_df %>%
  group_by(valid_id) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$valid_id))) %>%
  map(rebuild_web_objects)

# Calculate intermediate reference:
# Re_intdiag = stability of A + D_interaction
# where A is the off-diagonal trophic matrix and D_interaction is the
# interaction-derived diagonal only, without D_background.
intdiag_reference <- map_dfr(web_objects, function(obj) {
  J_intdiag <- build_jacobian(
    obj$C, obj$B, obj$e_link,
    J_diag = get_interaction_diag_ode(obj$C, obj$B, obj$e_link)
  )
  
  tibble(Re_intdiag = get_Re_lmax(J_intdiag))
}, .id = "valid_id")


web_meta <- modelfw_df %>%
  distinct(
    valid_id, S_label, rep,
    richness_eq, connectance_eq, maxTL_eq, L_eq,
    Re_emp, Re_ode
  ) %>%
  left_join(intdiag_reference, by = "valid_id") %>%
  mutate(
    reference_Re = reference_type,
    Re_reference = .data[[reference_type]]
  )

#------------------------------------------------------------
# 3. Validate rebuilt empirical, intdiag, and ODE-style Jacobians
#------------------------------------------------------------
validation_check <- map_dfr(names(web_objects), function(id) {
  obj <- web_objects[[id]]
  meta <- filter(web_meta, valid_id == id)
  
  J_emp <- build_jacobian(obj$C, obj$B, obj$e_link)
  
  J_intdiag <- build_jacobian(
    obj$C, obj$B, obj$e_link,
    J_diag = get_interaction_diag_ode(obj$C, obj$B, obj$e_link)
  )
  
  J_ode <- build_jacobian(
    obj$C, obj$B, obj$e_link,
    J_diag = obj$J_ode_diag
  )
  
  J_ode_off <- J_ode
  diag(J_ode_off) <- 0
  
  tibble(
    valid_id = id,
    diff_Re_emp = get_Re_lmax(J_emp) - meta$Re_emp,
    diff_Re_intdiag = get_Re_lmax(J_intdiag) - meta$Re_intdiag,
    diff_Re_ode = get_Re_lmax(J_ode) - meta$Re_ode,
    max_abs_J_emp_diff = max(abs(J_emp - obj$J_emp)),
    max_abs_offdiag_diff = max(abs(J_emp - J_ode_off)),
    min_J_ode_diag = min(obj$J_ode_diag),
    max_J_ode_diag = max(obj$J_ode_diag)
  )
})

max(abs(validation_check$diff_Re_emp))
max(abs(validation_check$diff_Re_intdiag))
max(validation_check$max_abs_J_emp_diff)
max(validation_check$max_abs_offdiag_diff)
max(abs(validation_check$diff_Re_ode))

summary(validation_check$diff_Re_emp)
summary(validation_check$diff_Re_intdiag)
summary(validation_check$diff_Re_ode)
summary(validation_check$max_abs_offdiag_diff)

#------------------------------------------------------------
# 4. Reference stability table
#------------------------------------------------------------
reference_results <- web_meta %>%
  select(
    valid_id, S_label, 
    richness_eq, connectance_eq, maxTL_eq, L_eq,
    reference_Re, Re_reference, Re_emp, Re_intdiag, Re_ode
  )
summary(reference_results$Re_ode)

#------------------------------------------------------------
# 5. Randomization
#------------------------------------------------------------
# Re_emp:
#   J_randomization = A_randomization, diagonal = 0
#
# Re_intdiag:
#   J_randomization = A_randomization + D_interaction_rand
#   background diagonal removed
#
# Re_ode:
#   J_randomization = A_randomization + D_interaction_randomization + D_background
# 
# randomisation_results <- map_dfr(names(web_objects), function(id) {
#   message("Running web: ", id)
#   
#   obj <- web_objects[[id]]
#   meta <- filter(web_meta, valid_id == id)
#   
#   background_diag <- switch(
#     reference_type,
#     Re_emp = NULL,
#     Re_intdiag = numeric(length(obj$B)),
#     Re_ode = obj$J_ode_diag -
#       get_interaction_diag_ode(obj$C, obj$B, obj$e_link)
#   )
#   
#   map_dfr(scenario_levels, function(scen) {
#     tibble(
#       rep = seq_len(n_rep),
#       valid_id = id,
#       S_label = meta$S_label,
#       S_target = meta$S_target,
#       richness_eq = meta$richness_eq,
#       connectance_eq = meta$connectance_eq,
#       maxTL_eq = meta$maxTL_eq,
#       L_eq = meta$L_eq,
#       reference_Re = reference_type,
#       scenario = scen,
#       Re_lmax = map_dbl(seq_len(n_rep), function(x) {
#         J_rand <- randomise_jacobian(
#           obj$C, obj$B, obj$e_link, scen,
#           background_diag = background_diag
#         )
#         get_Re_lmax(J_rand)
#       })
#     )
#   })
# }) %>%
#   mutate(
#     scenario = factor(scenario, levels = scenario_levels),
#     S_label = factor(S_label, levels = S_label_levels)
#   )
# 
# saveRDS(randomisation_results, randomisation_rds)

randomisation_results <- readRDS(randomisation_rds)

#------------------------------------------------------------
# 6. Calculate web-level SES deviations
#------------------------------------------------------------
rep_check <- randomisation_results %>%
  count(valid_id, scenario, name = "n_realisations")
summary(rep_check$n_realisations)

scenario_summary_by_web <- randomisation_results %>%
  group_by(valid_id, scenario) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    sd_Re = sd(Re_lmax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    web_meta %>%
      select(
        valid_id, S_label, richness_eq, connectance_eq, maxTL_eq,
        L_eq, Re_emp, Re_intdiag, Re_ode,
        reference_Re, Re_reference
      ),
    by = "valid_id"
  ) %>%
  mutate(
    ses_deviation_from_reference =
      (mean_Re - Re_reference) / sd_Re,
    scenario = factor(scenario, levels = scenario_levels),
    S_label = factor(S_label, levels = S_label_levels)
  )


#------------------------------------------------------------
# 7. stabilising effect of full dependency
#------------------------------------------------------------
# For panels B-D and the LM:
#   full_dependency_effect =
#     (mean_Re(no_dependency) - mean_Re(full_dependency)) /
#     sd_Re(no_dependency)
# Positive values mean that removing dependency increases Re(lambda_max)

full_dep_df <- scenario_summary_by_web %>%
  filter(scenario %in% c("full_dependency", "no_dependency")) %>%
  mutate(
    z_log_richness = make_z(log(richness_eq)), #using ln richness
    z_log_connectance = make_z(log(connectance_eq)), #using ln connectance
    z_log_complexity = make_z(log(L_eq / richness_eq)), #using ln complexit
    z_maxTL = make_z(maxTL_eq)
  ) %>%
  select(
    valid_id, scenario, mean_Re, sd_Re,
    richness_eq, connectance_eq, maxTL_eq,
    z_log_richness, z_log_connectance, z_maxTL
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(mean_Re, sd_Re)
  ) %>%
  mutate(
    full_dependency_effect_raw =
      mean_Re_no_dependency - mean_Re_full_dependency,
    
    full_dependency_effect_std =
      full_dependency_effect_raw / sd_Re_no_dependency
  ) 

full_dep_lm <- lm(
  full_dependency_effect_std ~ z_log_richness + z_log_connectance + z_maxTL,
  data = full_dep_df
)


summary(full_dep_lm)
anova(full_dep_lm)

effect_table <- tidy(full_dep_lm, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    predictor = recode(
      term,
      z_log_richness = "Log richness",
      z_log_connectance = "Log connectance",
      z_maxTL = "Maximum trophic level"
    )
  ) %>%
  select(predictor, estimate, conf.low, conf.high, p.value)

model_r2 <- tibble(
  r_squared = summary(full_dep_lm)$r.squared,
  adj_r_squared = summary(full_dep_lm)$adj.r.squared
)

effect_table
model_r2

#------------------------------------------------------------
# 8. Model-adjusted predictions
#------------------------------------------------------------
make_prediction_df <- function(model, data, focal_var, label) {
  pred_df <- tibble(
    z_log_richness = 0,
    z_log_connectance = 0,
    z_maxTL = 0
  )
  
  focal_seq <- seq(
    min(data[[focal_var]], na.rm = TRUE),
    max(data[[focal_var]], na.rm = TRUE),
    length.out = 100
  )
  
  pred_df <- pred_df[rep(1, length(focal_seq)), ]
  pred_df[[focal_var]] <- focal_seq
  
  pred <- predict(model, newdata = pred_df, se.fit = TRUE)
  
  pred_df %>%
    mutate(
      x = .data[[focal_var]],
      fit = pred$fit,
      se = pred$se.fit,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      complexity = label
    )
}

pred_richness <- make_prediction_df(
  full_dep_lm,
  full_dep_df,
  "z_log_richness",
  "Richness"
)

pred_connectance <- make_prediction_df(
  full_dep_lm,
  full_dep_df,
  "z_log_connectance",
  "Connectance"
)

pred_maxTL <- make_prediction_df(
  full_dep_lm,
  full_dep_df,
  "z_maxTL",
  "Max trophic level"
)

slope_labs <- tidy(full_dep_lm, conf.int = TRUE) %>%
  filter(term %in% c("z_log_richness", "z_log_connectance", "z_maxTL")) %>%
  mutate(
    label = paste0(
      "β = ", sprintf("%.2f", estimate),
      "\n95% CI [", sprintf("%.2f", conf.low), ", ",
      sprintf("%.2f", conf.high), "]"
    )
  )

full_dep_ylim <- range(
  c(
    full_dep_df$full_dependency_effect_std,
    pred_richness$lower, pred_richness$upper,
    pred_connectance$lower, pred_connectance$upper,
    pred_maxTL$lower, pred_maxTL$upper
  ),
  na.rm = TRUE
)

full_dep_ylim <- full_dep_ylim + c(-0.05, 0.05) * diff(full_dep_ylim)

#------------------------------------------------------------
# 9. Plots
#------------------------------------------------------------
theme_dep <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
      axis.line = element_blank(),
      axis.text = element_text(colour = "black", size = base_size - 2),
      axis.title = element_text(colour = "black", size = base_size),
      axis.ticks = element_line(linewidth = 0.45),
      legend.position = "bottom",
     # legend.title = element_blank(),
      legend.text = element_text(size = base_size - 4),
      plot.tag = element_text(face = "bold", size = base_size)
    )
}

# Mixed model for SES
mod_letters <- lmerTest::lmer(
  ses_deviation_from_reference ~ S_label * scenario + (1 | valid_id),
  data = scenario_summary_by_web,
  contrasts = list(
    S_label = "contr.sum",
    scenario = "contr.sum"
  )
)

# Tests of richness bin, scenario, and their interaction
anova_result <- anova(
  mod_letters,
  type = 3,
  ddf = "Satterthwaite"
)
print(anova_result)

# Tukey comparisons separately within each richness bin
emm_scen <- emmeans(
  mod_letters,
  ~ scenario | S_label,
  lmer.df = "satterthwaite"
)

pairwise_tests <- summary(
  pairs(emm_scen, adjust = "tukey"),
  infer = c(TRUE, TRUE)
)
print(pairwise_tests)

# Comparison letters
scenario_letters <- multcomp::cld(
  emm_scen,
  by = "S_label",
  adjust = "tukey",
  Letters = letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  mutate(.group = gsub(" ", "", .group)) %>%
  select(S_label, scenario, .group)

# Position letters above the points
y_gap <- 0.04 * max(
  1,
  diff(range(
    scenario_summary_by_web$ses_deviation_from_reference,
    na.rm = TRUE
  ))
)

scenario_letters <- scenario_summary_by_web %>%
  group_by(S_label, scenario) %>%
  summarise(
    y = max(ses_deviation_from_reference, na.rm = TRUE) + y_gap,
    .groups = "drop"
  ) %>%
  left_join(scenario_letters, by = c("S_label", "scenario"))


p_box <- ggplot(
  scenario_summary_by_web,
  aes(
    S_label, ses_deviation_from_reference,
    fill = scenario, colour = scenario)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45") +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 0.5,
    outlier.shape = NA ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      jitter.height = 0,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.5) +
  geom_text(
    data = scenario_letters,
    aes(S_label, y, label = .group, group = scenario),
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    colour = "black",
    size = 6,
    vjust = 0) +
  scale_fill_manual(
    values = scenario_cols,
    labels = scenario_labs ) +
  scale_colour_manual(
    values = scenario_cols,
    labels = scenario_labs) +
  scale_x_discrete(labels = S_label_labs) +
  scale_y_continuous(limits = c(-2.3, 2.5),
                     expand = expansion(mult = c(0.04, 0.10))) +
  # scale_y_continuous(limits = c(-1, 6.5),
  #                    expand = expansion(mult = c(0.04, 0.10))) +
  labs(x = "Species richness", y = "SES", fill = "") +
  guides(
    colour = "none",
    fill = guide_legend(nrow = 1) ) +
  theme_dep(20)

p_box




# bottom panel
make_partial_plot <- function(data, model, xvar, xlab, show_y = FALSE) {
  
  beta_df <- broom::tidy(model, conf.int = TRUE) %>%
    filter(term == xvar)
  
  beta <- beta_df$estimate
  
  panel_title <- slope_labs$label[slope_labs$term == xvar]
  
  partial_df <- data %>%
    transmute(
      x = .data[[xvar]],
      y = residuals(model) + beta * .data[[xvar]]
    )
  
  line_df <- tibble(
    x = seq(
      min(partial_df$x, na.rm = TRUE),
      max(partial_df$x, na.rm = TRUE),
      length.out = 100
    ),
    fit   = beta_df$estimate * x,
    lower1 = beta_df$conf.low * x,
    upper1 = beta_df$conf.high * x
  ) %>%
    mutate(
      lower = pmin(lower1, upper1),
      upper = pmax(lower1, upper1)
    )
  
  ggplot(partial_df, aes(x, y)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey45"
    ) +
    geom_point(
      size = 1.7,
      alpha = 0.7,
      colour = "grey35"
    ) +
    geom_ribbon(
      data = line_df,
      aes(x = x, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      alpha = 0.22
    ) +
    geom_line(
      data = line_df,
      aes(x = x, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    labs(
      title = panel_title,
      x = xlab,
      y = if (show_y) {
        "Partial effect on standardised\nstabilising effect"
      } else {
        NULL
      }
    ) +
    theme_dep(20) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0,
        margin = margin(b = 8)
      )
    )
}

p_richness <- make_partial_plot(
  full_dep_df, full_dep_lm,
  "z_log_richness",
  "Ln(Richness),\nstandardised",
  TRUE
)

p_connectance <- make_partial_plot(
  full_dep_df, full_dep_lm,
  "z_log_connectance",
  "Ln(connectance),\nstandardised"
)

p_maxTL <- make_partial_plot(
  full_dep_df, full_dep_lm,
  "z_maxTL",
  "Max. trophic level,\nstandardised"
)

main_plot <- (p_box / (p_richness | p_connectance | p_maxTL)) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

main_plot

ggsave(fig_jpeg, main_plot, width = 15, height = 12, dpi = 600)
ggsave(fig_pdf, main_plot, width = 15, height = 12)


#------------------------------------------------------------
# 10. Diagnostics: diagonal and off-diagonal elements
#------------------------------------------------------------
jacobian_diagnostics <- map_dfr(names(web_objects), function(id) {
  obj <- web_objects[[id]]
  
  J_ode <- build_jacobian(
    obj$C, obj$B, obj$e_link,
    J_diag = obj$J_ode_diag
  )
  
  diagonal <- diag(J_ode)
  off_diagonal <- J_ode[row(J_ode) != col(J_ode)]
  off_diagonal <- off_diagonal[off_diagonal != 0]
  
  tibble(
    valid_id = id,
    median_abs_diagonal = median(abs(diagonal), na.rm = TRUE),
    median_abs_off_diagonal = median(abs(off_diagonal), na.rm = TRUE),
    median_abs_positive_off_diagonal =
      median(abs(off_diagonal[off_diagonal > 0]), na.rm = TRUE),
    median_abs_negative_off_diagonal =
      median(abs(off_diagonal[off_diagonal < 0]), na.rm = TRUE)
  )
}) %>%
  left_join(
    web_meta %>%
      select(valid_id, S_label, richness_eq, connectance_eq, maxTL_eq, Re_ode),
    by = "valid_id"
  )

summary(jacobian_diagnostics)

p_diag_offdiag <- ggplot(
  jacobian_diagnostics,
  aes(median_abs_diagonal, median_abs_off_diagonal, colour = maxTL_eq)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = 2.5, alpha = 0.55) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_viridis_c(option = "plasma", direction = -1) +
  #scale_colour_viridis_d(labels = S_label_labs) +
  labs(
    x = "Median absolute diagonal element",
    y = "Median absolute \nnon-zero \noff-diagonal element",
    colour = "maxTL"
  ) +
  theme_dep(20)

p_diag_offdiag

ggsave(
  here("figures", paste0("model_fws_diag_offdiag_maxTL", reference_suffix, ".jpeg")),
  p_diag_offdiag, width = 8, height = 7, dpi = 600
)

#------------------------------------------------------------
# 10. Diagnostics: diagonal and off-diagonal elements
#------------------------------------------------------------

# Matrix convention:
# C[i, j] > 0 means species j consumes species i
#
# Therefore:
#   colSums(C != 0) == 0 : basal species
#   rowSums(C != 0) == 0 : top consumer
#   otherwise            : intermediate consumer


#------------------------------------------------------------
# 10.1 Individual diagonal elements
#------------------------------------------------------------

jacobian_diagonal_elements <- map_dfr(names(web_objects), function(id) {
  
  obj <- web_objects[[id]]
  
  J_ode <- build_jacobian(
    obj$C, obj$B, obj$e_link,
    J_diag = obj$J_ode_diag
  )
  
  diagonal <- diag(J_ode)
  
  row_links <- rowSums(obj$C != 0)  # number of predators
  col_links <- colSums(obj$C != 0)  # number of resources
  
  trophic_role <- case_when(
    col_links == 0 & row_links > 0 ~ "Basal",
    col_links > 0 & row_links > 0 ~ "Intermediate consumer",
    col_links > 0 & row_links == 0 ~ "Top consumer",
    col_links == 0 & row_links == 0 ~ "Isolated"
  )
  
  tibble(
    valid_id = id,
    species_id = seq_along(diagonal),
    trophic_role = trophic_role,
    row_links = row_links,
    col_links = col_links,
    Jii = diagonal,
    abs_Jii = abs(diagonal)
  )
}) %>%
  left_join(
    web_meta %>%
      select(
        valid_id, S_label, richness_eq,
        connectance_eq, maxTL_eq, Re_ode
      ),
    by = "valid_id"
  )


#------------------------------------------------------------
# 10.2 Summary of Jii by trophic role
#------------------------------------------------------------

Jii_role_summary <- jacobian_diagonal_elements %>%
  group_by(trophic_role) %>%
  summarise(
    n = n(),
    median_abs_Jii = median(abs_Jii, na.rm = TRUE),
    q25 = quantile(abs_Jii, 0.25, na.rm = TRUE),
    q75 = quantile(abs_Jii, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

Jii_role_summary


#------------------------------------------------------------
# 10.3 Distribution of individual Jii
#------------------------------------------------------------

p_Jii_role <- jacobian_diagonal_elements %>%
  filter(
    is.finite(abs_Jii),
    abs_Jii > 0,
    trophic_role != "Isolated"
  ) %>%
  ggplot(
    aes(
      x = abs_Jii,
      fill = trophic_role
    )
  ) +
  geom_histogram(
    bins = 80,
    position = "identity",
    alpha = 0.55
  ) +
  scale_x_log10() +
  labs(
    x = expression(abs(J[ii])),
    y = "Number of species",
    fill = "Trophic role"
  ) +
  theme_dep(20)

p_Jii_role
ggsave(
  here("figures", paste0("model_fws_Jii_distribution", reference_suffix, ".jpeg")),
  p_Jii_role, width = 10, height = 7, dpi = 600
)

#------------------------------------------------------------
# 10.4 Food-web trophic composition
#------------------------------------------------------------

web_trophic_composition <- jacobian_diagonal_elements %>%
  group_by(valid_id) %>%
  summarise(
    n_species = n(),
    
    n_basal =
      sum(trophic_role == "Basal"),
    
    n_intermediate =
      sum(trophic_role == "Intermediate consumer"),
    
    n_top =
      sum(trophic_role == "Top consumer"),
    
    prop_basal =
      mean(trophic_role == "Basal"),
    
    prop_intermediate =
      mean(trophic_role == "Intermediate consumer"),
    
    prop_top =
      mean(trophic_role == "Top consumer"),
    
    .groups = "drop"
  )


#------------------------------------------------------------
# 10.5 Web-level diagonal and off-diagonal summaries
#------------------------------------------------------------

jacobian_diagnostics <- map_dfr(names(web_objects), function(id) {
  
  obj <- web_objects[[id]]
  
  J_ode <- build_jacobian(
    obj$C, obj$B, obj$e_link,
    J_diag = obj$J_ode_diag
  )
  
  diagonal <- diag(J_ode)
  
  off_diagonal <- J_ode[
    row(J_ode) != col(J_ode)
  ]
  
  off_diagonal <- off_diagonal[
    off_diagonal != 0
  ]
  
  tibble(
    valid_id = id,
    
    median_abs_diagonal =
      median(abs(diagonal), na.rm = TRUE),
    
    median_abs_off_diagonal =
      median(abs(off_diagonal), na.rm = TRUE),
    
    median_abs_positive_off_diagonal =
      median(
        abs(off_diagonal[off_diagonal > 0]),
        na.rm = TRUE
      ),
    
    median_abs_negative_off_diagonal =
      median(
        abs(off_diagonal[off_diagonal < 0]),
        na.rm = TRUE
      )
  )
}) %>%
  
  left_join(
    web_meta %>%
      select(
        valid_id, S_label, richness_eq,
        connectance_eq, maxTL_eq, Re_ode
      ),
    by = "valid_id"
  ) %>%
  
  left_join(
    web_trophic_composition,
    by = "valid_id"
  )

summary(jacobian_diagnostics)


#------------------------------------------------------------
# 10.6 Median diagonal vs off-diagonal
#      coloured by proportion of basal species
#------------------------------------------------------------
basal_range <- range(
  jacobian_diagnostics$prop_basal,
  na.rm = TRUE
)

p_diag_offdiag <- ggplot(
  jacobian_diagnostics,
  aes(
    median_abs_diagonal,
    median_abs_off_diagonal,
    colour = prop_basal
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey45"
  ) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_viridis_c(
    option = "cividis",
    limits = basal_range,
    guide = guide_colourbar(
      barwidth = unit(5.5, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(
    x = "Median absolute diagonal element",
    y = "Median absolute non-zero\noff-diagonal element",
    colour = "Proportion\nbasal"
  ) +
  theme_dep(20)

p_diag_offdiag

ggsave(
  here("figures", paste0("model_fws_diag_offdiag", reference_suffix, ".jpeg")),
  p_diag_offdiag, width = 8, height = 7, dpi = 600
)
