#============================================================
# Whole-network dependency randomisation for niche-model food webs
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

select <- dplyr::select
filter <- dplyr::filter

set.seed(2026)

source(here("R/internal.R"))

#------------------------------------------------------------
# 1. Settings
#------------------------------------------------------------

out_dir <- here("data/model_fws/linear_response_sim")

feather_file <- file.path(out_dir, "model_food_webs.feather")

# Choose baseline stability:
# "Re_emp"     = empirical-style Jacobian, diagonal = 0
# "Re_intdiag" = off-diagonal trophic terms + interaction-derived diagonal only
# "Re_ode"     = full ODE Jacobian, including interaction-derived and background diagonal

# reference_type <- "Re_emp"
 reference_type <- "Re_intdiag"
# reference_type <- "Re_ode"

reference_label <- recode(
  reference_type,
  Re_emp     = "Re(emp)",
  Re_intdiag = "Re(int. diag)",
  Re_ode     = "Re(ode)"
)

reference_suffix <- paste0("_", reference_type)

randomisation_rds <- file.path(
  out_dir,
  paste0("model_fw_dependency_randomisation", reference_suffix, ".rds")
)

fig_jpeg <- here("figures", paste0("model_fws_dependency_combined", reference_suffix, ".jpeg"))
fig_pdf  <- here("figures", paste0("model_fws_dependency_combined", reference_suffix, ".pdf"))

n_rep <- 1000

S_label_levels <- c("6_10", "11_15", "16_20", "21_25", "26_30")

S_label_labs <- c(
  "6_10"  = "6–10",
  "11_15" = "11–15",
  "16_20" = "16–20",
  "21_25" = "21–25",
  "26_30" = "26–30"
)

scenario_levels <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

scenario_cols <- c(
  full_dependency = "#1B9E77",
  c_dependency    = "#7570B3",
  B_dependency    = "#E6AB02",
  e_dependency    = "#D95F02",
  no_dependency   = "#666666"
)

scenario_cols_fill <- scales::alpha(scenario_cols, 0.25)

scenario_labs <- c(
  full_dependency = "Full structural\ndependency",
  c_dependency    = "Interaction coefficient\nc dependency",
  B_dependency    = "Shared biomass\nB dependency",
  e_dependency    = "Shared efficiency\ne dependency",
  no_dependency   = "No dependency"
)

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
intdiag_reference <- map_dfr(names(web_objects), function(id) {
  obj <- web_objects[[id]]
  
  interaction_diag_original <- get_interaction_diag_ode(
    obj$C,
    obj$B,
    obj$e_link
  )
  
  J_intdiag_original <- build_jacobian(
    obj$C,
    obj$B,
    obj$e_link,
    J_diag = interaction_diag_original
  )
  
  tibble(
    valid_id = id,
    Re_intdiag = get_Re_lmax(J_intdiag_original)
  )
})

web_meta <- modelfw_df %>%
  distinct(
    valid_id, S_label, S_target, rep,
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
  meta <- web_meta %>% filter(valid_id == id)
  
  Re_emp_ref <- unique(meta$Re_emp)
  Re_intdiag_ref <- unique(meta$Re_intdiag)
  Re_ode_ref <- unique(meta$Re_ode)
  
  stopifnot(length(Re_emp_ref) == 1)
  stopifnot(length(Re_intdiag_ref) == 1)
  stopifnot(length(Re_ode_ref) == 1)
  
  # Empirical-style Jacobian: A only, diagonal = 0
  J_emp_rebuilt <- build_jacobian(obj$C, obj$B, obj$e_link)
  
  # Interaction-derived diagonal only: A + D_interaction
  interaction_diag_original <- get_interaction_diag_ode(
    obj$C,
    obj$B,
    obj$e_link
  )
  
  J_intdiag_rebuilt <- build_jacobian(
    obj$C,
    obj$B,
    obj$e_link,
    J_diag = interaction_diag_original
  )
  
  diff_Re_ode <- NA_real_
  min_J_ode_diag <- NA_real_
  max_J_ode_diag <- NA_real_
  max_abs_offdiag_diff <- NA_real_
  
  if (!is.null(obj$J_ode_diag)) {
    J_ode_rebuilt <- build_jacobian(
      obj$C,
      obj$B,
      obj$e_link,
      J_diag = obj$J_ode_diag
    )
    
    J_ode_off <- J_ode_rebuilt
    diag(J_ode_off) <- 0
    
    max_abs_offdiag_diff <- max(abs(J_emp_rebuilt - J_ode_off), na.rm = TRUE)
    
    diff_Re_ode <- get_Re_lmax(J_ode_rebuilt) - as.numeric(Re_ode_ref)
    min_J_ode_diag <- min(obj$J_ode_diag, na.rm = TRUE)
    max_J_ode_diag <- max(obj$J_ode_diag, na.rm = TRUE)
  }
  
  tibble(
    valid_id = id,
    diff_Re_emp = get_Re_lmax(J_emp_rebuilt) - as.numeric(Re_emp_ref),
    diff_Re_intdiag = get_Re_lmax(J_intdiag_rebuilt) - as.numeric(Re_intdiag_ref),
    max_abs_J_emp_diff = max(abs(J_emp_rebuilt - obj$J_emp), na.rm = TRUE),
    max_abs_offdiag_diff = max_abs_offdiag_diff,
    diff_Re_ode = diff_Re_ode,
    min_J_ode_diag = min_J_ode_diag,
    max_J_ode_diag = max_J_ode_diag
  )
})

print(max(abs(validation_check$diff_Re_emp), na.rm = TRUE))
print(max(abs(validation_check$diff_Re_intdiag), na.rm = TRUE))
print(max(validation_check$max_abs_J_emp_diff, na.rm = TRUE))
print(max(validation_check$max_abs_offdiag_diff, na.rm = TRUE))
print(max(abs(validation_check$diff_Re_ode), na.rm = TRUE))

summary(validation_check$diff_Re_emp)
summary(validation_check$diff_Re_intdiag)
summary(validation_check$diff_Re_ode)
summary(validation_check$max_abs_offdiag_diff)

#------------------------------------------------------------
# 4. Reference stability table
#------------------------------------------------------------

reference_results <- web_meta %>%
  select(
    valid_id, S_label, S_target,
    richness_eq, connectance_eq, maxTL_eq, L_eq,
    reference_Re, Re_reference, Re_emp, Re_intdiag, Re_ode
  )

#------------------------------------------------------------
# 5. Randomisation
#------------------------------------------------------------
# Important:
# Re_emp:
#   J_rand = A_rand, diagonal = 0
#
# Re_intdiag:
#   J_rand = A_rand + D_interaction_rand
#   background diagonal removed
#
# Re_ode:
#   J_rand = A_rand + D_interaction_rand + D_background

randomisation_results <- map_dfr(names(web_objects), function(id) {
  message("Running web: ", id)
  
  obj <- web_objects[[id]]
  meta <- web_meta %>% filter(valid_id == id)
  
  if (identical(reference_type, "Re_emp")) {
    
    background_diag <- NULL
    
  } else if (identical(reference_type, "Re_intdiag")) {
    
    if (is.null(obj$J_ode_diag)) {
      stop("J_ode_diag is missing. Update the Feather file before running Re_intdiag.")
    }
    
    # Remove D_background.
    # compute_randomised_Re_lmax_ode() will still add D_interaction_rand.
    background_diag <- rep(0, length(obj$J_ode_diag))
    
  } else if (identical(reference_type, "Re_ode")) {
    
    if (is.null(obj$J_ode_diag)) {
      stop("J_ode_diag is missing. Update the Feather file before running Re_ode.")
    }
    
    interaction_diag_original <- get_interaction_diag_ode(
      obj$C,
      obj$B,
      obj$e_link
    )
    
    background_diag <- obj$J_ode_diag - interaction_diag_original
    
  } else {
    stop("Unknown reference_type: ", reference_type)
  }
  
  map_dfr(scenario_levels, function(scen) {
    tibble(
      rep = seq_len(n_rep),
      valid_id = id,
      S_label = meta$S_label,
      S_target = meta$S_target,
      richness_eq = meta$richness_eq,
      connectance_eq = meta$connectance_eq,
      maxTL_eq = meta$maxTL_eq,
      L_eq = meta$L_eq,
      reference_Re = reference_type,
      scenario = scen,
      Re_lmax = map_dbl(
        seq_len(n_rep),
        ~ if (identical(reference_type, "Re_emp")) {
          compute_randomised_Re_lmax(
            obj$C,
            obj$B,
            obj$e_link,
            scen
          )
        } else {
          compute_randomised_Re_lmax_ode(
            obj$C,
            obj$B,
            obj$e_link,
            scen,
            background_diag
          )
        }
      )
    )
  })
}) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    S_label = factor(S_label, levels = S_label_levels)
  )

saveRDS(randomisation_results, randomisation_rds)

randomisation_results <- readRDS(randomisation_rds)

#------------------------------------------------------------
# 6. Calculate web-level SES deviations
#------------------------------------------------------------

rep_check <- randomisation_results %>%
  count(valid_id, scenario, name = "n_realisations")

print(summary(rep_check$n_realisations))

rep_check %>%
  count(n_realisations) %>%
  print()

rep_check %>%
  filter(n_realisations != n_rep) %>%
  print()

scenario_summary_by_web <- randomisation_results %>%
  group_by(
    valid_id, S_label, S_target,
    richness_eq, connectance_eq, maxTL_eq, L_eq,
    scenario
  ) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    sd_Re = safe_sd(Re_lmax),
    .groups = "drop"
  ) %>%
  left_join(
    web_meta %>%
      dplyr::select(valid_id, reference_Re, Re_reference, Re_emp, Re_intdiag, Re_ode),
    by = "valid_id"
  )

analysis_df <- scenario_summary_by_web %>%
  mutate(
    raw_deviation_from_reference = mean_Re - Re_reference,
    
    # Panel A:
    # scenario-specific SES deviation from reference stability
    ses_deviation_from_reference =
      raw_deviation_from_reference / sd_Re,
    
    z_log_richness = make_z(log(richness_eq)),
    z_log_connectance = make_z(log(connectance_eq)),
    z_maxTL = make_z(maxTL_eq),
    
    scenario = factor(scenario, levels = scenario_levels),
    S_label = factor(S_label, levels = S_label_levels)
  ) %>%
  filter(
    is.finite(ses_deviation_from_reference),
    is.finite(z_log_richness),
    is.finite(z_log_connectance),
    is.finite(z_maxTL)
  )

# Optional check for very small SDs
scenario_summary_by_web %>%
  group_by(scenario) %>%
  summarise(
    min_sd_Re = min(sd_Re, na.rm = TRUE),
    median_sd_Re = median(sd_Re, na.rm = TRUE),
    n_zero_sd = sum(sd_Re == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

#------------------------------------------------------------
# 7. Main analysis: stabilising effect of full dependency
#------------------------------------------------------------
# For panels B-D and the LM:
#
#   full_dependency_effect =
#     (mean_Re(no_dependency) - mean_Re(full_dependency)) /
#     sd_Re(no_dependency)
#
# Positive values mean that removing dependency increases Re(lambda_max),
# so full dependency is stabilising.

full_dep_df <- scenario_summary_by_web %>%
  filter(scenario %in% c("full_dependency", "no_dependency")) %>%
  mutate(
    z_log_richness = make_z(log(richness_eq)),
    z_log_connectance = make_z(log(connectance_eq)),
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
    
    full_dependency_effect =
      full_dependency_effect_raw / sd_Re_no_dependency
  ) %>%
  filter(
    is.finite(full_dependency_effect),
    is.finite(z_log_richness),
    is.finite(z_log_connectance),
    is.finite(z_maxTL)
  )

full_dep_lm <- lm(
  full_dependency_effect ~ z_log_richness + z_log_connectance + z_maxTL,
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

print(effect_table)
print(model_r2)

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
    full_dep_df$full_dependency_effect,
    pred_richness$lower, pred_richness$upper,
    pred_connectance$lower, pred_connectance$upper,
    pred_maxTL$lower, pred_maxTL$upper
  ),
  na.rm = TRUE
)

full_dep_ylim <- full_dep_ylim + c(-0.05, 0.05) * diff(full_dep_ylim)

#------------------------------------------------------------
# 9. Combined figure
#------------------------------------------------------------
theme_dep <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black", size = base_size - 2),
      axis.title = element_text(colour = "black", size = base_size),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.45),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2),
      plot.tag = element_text(face = "bold", size = base_size, colour = "black")
    )
}

plot_df <- analysis_df %>%
  filter(!is.na(S_label), is.finite(ses_deviation_from_reference)) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    S_label = factor(S_label, levels = S_label_levels)
  )

#------------------------------------------------------------
# Pairwise comparisons among scenarios within each richness bin
#------------------------------------------------------------
mod_letters <- lmer(
  ses_deviation_from_reference ~ S_label * scenario + (1 | valid_id),
  data = plot_df
)

emm_scen <- emmeans(
  mod_letters,
  specs = ~ scenario,
  by = "S_label"
)

pairwise_scen <- contrast(
  emm_scen,
  method = "pairwise",
  by = "S_label",
  adjust = "tukey"
)

print(pairwise_scen)

scenario_letters <- multcomp::cld(
  emm_scen,
  by = "S_label",
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  mutate(
    .group = gsub(" ", "", .group),
    scenario = factor(scenario, levels = scenario_levels),
    S_label = factor(S_label, levels = S_label_levels)
  )

letter_pos <- plot_df %>%
  group_by(S_label, scenario) %>%
  summarise(
    y = quantile(ses_deviation_from_reference, 0.98, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(S_label) %>%
  mutate(
    y = y + 0.25 * diff(range(plot_df$ses_deviation_from_reference, na.rm = TRUE))
  ) %>%
  ungroup()

scenario_letters <- scenario_letters %>%
  left_join(letter_pos, by = c("S_label", "scenario"))

scenario_letters %>%
  dplyr::select(
    dplyr::any_of(c(
      "S_label", "scenario", "emmean", "SE", "df",
      "lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL", ".group"
    ))
  ) %>%
  tibble::as_tibble() %>%
  print(n = Inf)

p_box <- ggplot(
  plot_df,
  aes(S_label, ses_deviation_from_reference, fill = scenario, colour = scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_boxplot(
    position = position_dodge(0.82),
    width = 0.58,
    alpha = 0.58,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.82),
    alpha = 1,
    size = 1.5,
    stroke = 0
  ) +
  geom_text(
    data = scenario_letters,
    aes(
      x = S_label,
      y = y,
      label = .group,
      group = scenario
    ),
    position = position_dodge(width = 0.82),
    inherit.aes = FALSE,
    colour = "black",
    size = 5.5,
    vjust = 0
  ) +
  scale_fill_manual(values = scenario_cols_fill, labels = scenario_labs) +
  scale_colour_manual(values = scenario_cols, labels = scenario_labs) +
  scale_x_discrete(labels = S_label_labs) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000),
    labels = scales::label_number()
  ) +
  labs(
    x = "Richness bin",
    y = paste0(
      "Scenario-specific SES deviation from\n",
      reference_label,
      " stability"
    )
  ) +
  guides(colour = "none", fill = guide_legend(nrow = 1)) +
  theme_dep(base_size = 20)

p_box

plot_complexity <- function(data, pred, xvar, xlab, term, show_y = FALSE) {
  lab_df <- slope_labs %>%
    filter(term == !!term) %>%
    mutate(
      x = min(data[[xvar]], na.rm = TRUE),
      y = full_dep_ylim[2] + 0.2 * diff(full_dep_ylim)
    )
  
  ggplot(data, aes(.data[[xvar]], full_dependency_effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_point(size = 1.5, alpha = 0.7, colour = "grey35") +
    geom_ribbon(
      data = pred,
      aes(x = x, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      alpha = 0.22
    ) +
    geom_line(
      data = pred,
      aes(x = x, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_text(
      data = lab_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 5.5
    ) +
    coord_cartesian(
      ylim = full_dep_ylim,
      clip = "off"
    ) +
    labs(
      x = xlab,
      y = ifelse(
        show_y,
        "Standardised stabilising effect\nof full dependency",
        ""
      )
    ) +
    theme_dep(base_size = 20) +
    theme(legend.position = "none")
}

p_richness <- plot_complexity(
  full_dep_df,
  pred_richness,
  "z_log_richness",
  "Log richness,\nstandardised",
  "z_log_richness",
  TRUE
)

p_connectance <- plot_complexity(
  full_dep_df,
  pred_connectance,
  "z_log_connectance",
  "Log connectance,\nstandardised",
  "z_log_connectance"
)

p_maxTL <- plot_complexity(
  full_dep_df,
  pred_maxTL,
  "z_maxTL",
  "Max. trophic level,\nstandardised",
  "z_maxTL"
)

main_plot <- (p_box / (p_richness | p_connectance | p_maxTL)) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 20))

main_plot

ggsave(fig_jpeg, main_plot, width = 14, height = 11, dpi = 600)
ggsave(fig_pdf, main_plot, width = 14, height = 11, device = cairo_pdf)


#------------------------------------------------------------
# 10. Alternative combined figure: partial residual plots for B-D
#------------------------------------------------------------
fig_jpeg_partial <- here("figures", paste0("model_fws_dependency_combined_partial", reference_suffix, ".jpeg"))
fig_pdf_partial  <- here("figures", paste0("model_fws_dependency_combined_partial", reference_suffix, ".pdf"))

make_partial_plot <- function(data, model, xvar, xlab, term, show_y = FALSE) {
  beta_df <- tidy(model, conf.int = TRUE) %>% filter(term == !!term)
  beta <- beta_df$estimate
  
  plot_df <- data %>%
    mutate(
      x = .data[[xvar]],
      partial_effect = residuals(model) + beta * .data[[xvar]]
    )
  
  line_df <- tibble(
    x = seq(min(plot_df$x, na.rm = TRUE), max(plot_df$x, na.rm = TRUE), length.out = 100)
  ) %>%
    mutate(
      fit = beta_df$estimate * x,
      lower_raw = beta_df$conf.low * x,
      upper_raw = beta_df$conf.high * x,
      lower = pmin(lower_raw, upper_raw),
      upper = pmax(lower_raw, upper_raw)
    )
  
  y_range <- range(
    c(plot_df$partial_effect, line_df$lower, line_df$upper),
    na.rm = TRUE
  )
  
  y_pad <- 0.2 * diff(y_range)
  y_lim <- y_range + c(-0.04, 0.04) * diff(y_range)
  
  lab_df <- slope_labs %>%
    filter(term == !!term) %>%
    mutate(
      x = min(plot_df$x, na.rm = TRUE),
      y = y_lim[2] + y_pad
    )
  
  ggplot(plot_df, aes(x, partial_effect)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_point(size = 1.7, alpha = 0.7, colour = "grey35") +
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
    geom_text(
      data = lab_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 5
    ) +
    coord_cartesian(
      ylim = y_lim,
      clip = "off"
    ) +
    labs(
      x = xlab,
      y = ifelse(
        show_y,
        "Partial effect on standardised\nstabilising effect",
        ""
      )
    ) +
    theme_dep(base_size = 20) +
    theme(
      legend.position = "none",
      plot.margin = margin(15, 10, 10, 10)
    )
}

p_richness_partial <- make_partial_plot(
  full_dep_df, full_dep_lm,
  "z_log_richness",
  "Log richness,\nstandardised",
  "z_log_richness",
  TRUE
)

p_connectance_partial <- make_partial_plot(
  full_dep_df, full_dep_lm,
  "z_log_connectance",
  "Log connectance,\nstandardised",
  "z_log_connectance"
)

p_maxTL_partial <- make_partial_plot(
  full_dep_df, full_dep_lm,
  "z_maxTL",
  "Max. trophic level,\nstandardised",
  "z_maxTL"
)

main_plot_partial <- (p_box / (p_richness_partial | p_connectance_partial | p_maxTL_partial)) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 20))

main_plot_partial

ggsave(fig_jpeg_partial, main_plot_partial, width = 14, height = 11, dpi = 600)
ggsave(fig_pdf_partial, main_plot_partial, width = 14, height = 11, device = cairo_pdf)

