#============================================================
# 3. Whole-network dependency randomization for niche modeled food webs
#============================================================
rm(list = ls())

library(tidyverse)
library(here)
library(arrow)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)
library(scico)
library(colorspace)
library(patchwork)

set.seed(2026)

#------------------------------------------------------------
# 1. Settings
#------------------------------------------------------------

n_rep <- 1000

scenarios <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

scenario_levels <- c(
  "full_dependency",
  "c_dependency",
  "B_dependency",
  "e_dependency",
  "no_dependency"
)

scenario_labs <- c(
  full_dependency = "Full structural\ndependency",
  c_dependency    = "Interaction coefficient\nc dependency",
  B_dependency    = "Shared biomass\nB dependency",
  e_dependency    = "Shared efficiency\ne dependency",
  no_dependency   = "No dependency"
)

#------------------------------------------------------------
# 2. Read Feather file from Julia
#------------------------------------------------------------

dep_df <- read_feather(
  here(
    "data/model_fws/linear_response_niche_webs",
    "valid_food_webs_dependency_data.feather"
  )
)

dep_df <- dep_df %>%
  mutate(
    valid_id = as.character(valid_id),
    S_label = factor(S_label, levels = c("low", "medium", "high"))
  )

# quick check
glimpse(dep_df)

#------------------------------------------------------------
# 3. Basic helper functions
#------------------------------------------------------------

get_Re_lmax <- function(J) {
  ev <- eigen(J, only.values = TRUE)$values
  max(Re(ev))
}

rebuild_web_objects <- function(df) {
  richness_eq <- unique(df$richness_eq)
  stopifnot(length(richness_eq) == 1)
  
  n <- richness_eq
  
  C <- matrix(0, n, n)
  J_emp <- matrix(0, n, n)
  
  B <- rep(NA_real_, n)
  e_link <- matrix(0, n, n)
  
  # Rows in C are prey/resource, columns are predator/consumer
  C[cbind(df$prey, df$predator)] <- df$cij
  
  # Rebuild biomass vector
  B[df$prey] <- df$B_prey
  B[df$predator] <- df$B_predator
  
  # Rebuild link-specific efficiency matrix
  # Rows are prey/resource, columns are predator/consumer
  e_link[cbind(df$prey, df$predator)] <- df$e
  
  # Rebuild J_emp from saved top-down and bottom-up terms
  J_emp[cbind(df$prey, df$predator)] <- df$J_topdown
  J_emp[cbind(df$predator, df$prey)] <- df$J_bottomup
  
  diag(J_emp) <- 0
  
  list(
    C = C,
    B = B,
    e_link = e_link,
    J_emp = J_emp
  )
}

build_jacobian <- function(C, B, e_link) {
  n <- length(B)
  J <- matrix(0, n, n)
  
  links <- which(C > 0, arr.ind = TRUE)
  
  if (nrow(links) > 0) {
    for (k in seq_len(nrow(links))) {
      i <- links[k, 1]  # prey/resource
      j <- links[k, 2]  # predator/consumer
      
      J[i, j] <- J[i, j] - C[i, j] * B[i]
      J[j, i] <- J[j, i] + e_link[i, j] * C[i, j] * B[j]
    }
  }
  
  diag(J) <- 0
  J
}

#------------------------------------------------------------
# 4. Rebuild original food-web objects
#------------------------------------------------------------

web_objects <- dep_df %>%
  group_by(valid_id) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$valid_id))) %>%
  map(rebuild_web_objects)

web_meta <- dep_df %>%
  distinct(
    valid_id,
    S_label,
    S_min,
    S_max,
    S_target,
    rep,
    richness_eq,
    connectance_eq,
    maxTL_eq,
    L_eq,
    Re_emp,
    Re_ode
  )

#------------------------------------------------------------
# 5. Validate rebuilt Jacobian against Julia J_emp
#------------------------------------------------------------

validation_results <- map_dfr(names(web_objects), function(id) {
  obj <- web_objects[[id]]
  
  J_rebuilt <- build_jacobian(
    C = obj$C,
    B = obj$B,
    e_link = obj$e_link
  )
  
  Re_rebuilt <- get_Re_lmax(J_rebuilt)
  Re_saved <- web_meta$Re_emp[match(id, web_meta$valid_id)]
  
  tibble(
    valid_id = id,
    Re_saved = Re_saved,
    Re_rebuilt = Re_rebuilt,
    diff_Re = Re_rebuilt - Re_saved,
    max_abs_J_diff = max(abs(J_rebuilt - obj$J_emp), na.rm = TRUE)
  )
})

summary(validation_results$diff_Re)
summary(validation_results$max_abs_J_diff)

validation_results

#------------------------------------------------------------
# 6. Parameter pools
#------------------------------------------------------------

get_parameter_pools <- function(C, B, e_link) {
  list(
    c_pool = as.numeric(C[C > 0]),
    B_pool = as.numeric(B[B > 0]),
    e_pool = as.numeric(e_link[e_link > 0])
  )
}

get_links <- function(C) {
  links <- which(C > 0, arr.ind = TRUE)
  
  tibble(
    prey = links[, 1],
    predator = links[, 2]
  )
}

#------------------------------------------------------------
# 7. Dependency randomisation
#------------------------------------------------------------

randomise_params <- function(C, B, e_link, scenario) {
  n <- length(B)
  
  pools <- get_parameter_pools(C, B, e_link)
  links <- get_links(C)
  
  c_pool <- pools$c_pool
  B_pool <- pools$B_pool
  e_pool <- pools$e_pool
  
  if (length(c_pool) == 0) stop("c_pool is empty.")
  if (length(B_pool) == 0) stop("B_pool is empty.")
  if (length(e_pool) == 0) stop("e_pool is empty.")
  
  #----------------------------------------------------------
  # Full dependency
  # One c per trophic link, one B per species,
  # and one e per realised trophic link.
  #----------------------------------------------------------
  if (scenario == "full_dependency") {
    C_new <- matrix(0, n, n)
    e_new <- matrix(0, n, n)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      C_new[i, j] <- sample(c_pool, 1, replace = TRUE)
      e_new[i, j] <- sample(e_pool, 1, replace = TRUE)
    }
    
    B_new <- sample(B_pool, n, replace = TRUE)
    
    return(list(C = C_new, B = B_new, e_link = e_new))
  }
  
  #----------------------------------------------------------
  # c dependency only
  # Same c_ij is reused in top-down and bottom-up terms.
  # B and e are independently sampled.
  #----------------------------------------------------------
  if (scenario == "c_dependency") {
    J <- matrix(0, n, n)
    
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
  # B dependency only
  # One biomass per species is reused.
  # c and e are independently sampled for each occurrence.
  #----------------------------------------------------------
  if (scenario == "B_dependency") {
    B_new <- sample(B_pool, n, replace = TRUE)
    J <- matrix(0, n, n)
    
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
  # e dependency only
  # One e value is reused per predator.
  # This is closest to your Lake Constance logic,
  # where one efficiency belongs to a consumer.
  #----------------------------------------------------------
  if (scenario == "e_dependency") {
    predators <- sort(unique(links$predator))
    
    e_by_predator <- rep(0, n)
    e_by_predator[predators] <- sample(
      e_pool,
      length(predators),
      replace = TRUE
    )
    
    J <- matrix(0, n, n)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      c_topdown <- sample(c_pool, 1, replace = TRUE)
      B_topdown <- sample(B_pool, 1, replace = TRUE)
      
      c_bottomup <- sample(c_pool, 1, replace = TRUE)
      B_bottomup <- sample(B_pool, 1, replace = TRUE)
      
      J[i, j] <- J[i, j] - c_topdown * B_topdown
      J[j, i] <- J[j, i] + e_by_predator[j] * c_bottomup * B_bottomup
    }
    
    diag(J) <- 0
    return(list(J = J))
  }
  
  #----------------------------------------------------------
  # No dependency
  # Every occurrence of c, B, and e is independently sampled.
  #----------------------------------------------------------
  if (scenario == "no_dependency") {
    J <- matrix(0, n, n)
    
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

compute_randomised_Re_lmax <- function(C, B, e_link, scenario) {
  rand <- randomise_params(C, B, e_link, scenario)
  
  if (!is.null(rand$J)) {
    return(get_Re_lmax(rand$J))
  }
  
  J <- build_jacobian(
    C = rand$C,
    B = rand$B,
    e_link = rand$e_link
  )
  
  get_Re_lmax(J)
}

#------------------------------------------------------------
# 8. Original results from saved J_emp
#------------------------------------------------------------

original_results <- map_dfr(names(web_objects), function(id) {
  obj <- web_objects[[id]]
  meta <- web_meta %>% filter(valid_id == id)
  
  tibble(
    valid_id = id,
    S_label = meta$S_label,
    S_target = meta$S_target,
    richness_eq = meta$richness_eq,
    connectance_eq = meta$connectance_eq,
    maxTL_eq = meta$maxTL_eq,
    L_eq = meta$L_eq,
    Re_original = get_Re_lmax(obj$J_emp),
    Re_emp_saved = meta$Re_emp,
    Re_ode = meta$Re_ode
  )
})

#------------------------------------------------------------
# 9. Randomisation loop
#------------------------------------------------------------

randomisation_results <- map_dfr(names(web_objects), function(id) {
  message("Running web: ", id)
  
  obj <- web_objects[[id]]
  meta <- web_meta %>% filter(valid_id == id)
  
  map_dfr(scenarios, function(scen) {
    message("  Scenario: ", scen)
    
    tibble(rep = seq_len(n_rep)) %>%
      mutate(
        valid_id = id,
        S_label = meta$S_label,
        S_target = meta$S_target,
        richness_eq = meta$richness_eq,
        connectance_eq = meta$connectance_eq,
        maxTL_eq = meta$maxTL_eq,
        L_eq = meta$L_eq,
        scenario = scen,
        Re_lmax = map_dbl(
          rep,
          ~ compute_randomised_Re_lmax(
            C = obj$C,
            B = obj$B,
            e_link = obj$e_link,
            scenario = scen
          )
        )
      )
  })
})

saveRDS(
  randomisation_results,
  here("data/model_fws/linear_response_niche_webs",
       "niche_dependency_randomisation_results.rds")
)

#------------------------------------------------------------
# 10. Summarise randomisation results
#------------------------------------------------------------

scenario_summary_by_web <- randomisation_results %>%
  group_by(
    valid_id,
    S_label,
    S_target,
    richness_eq,
    connectance_eq,
    maxTL_eq,
    L_eq,
    scenario
  ) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    sd_Re = sd(Re_lmax, na.rm = TRUE),
    median_Re = median(Re_lmax, na.rm = TRUE),
    q025_Re = quantile(Re_lmax, 0.025, na.rm = TRUE),
    q975_Re = quantile(Re_lmax, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels)
  )

effect_size_results <- scenario_summary_by_web %>%
  left_join(
    original_results %>%
      select(valid_id, Re_original),
    by = "valid_id"
  ) %>%
  mutate(
    raw_effect = mean_Re - Re_original,
    standardised_effect = (Re_original - mean_Re) / sd_Re
  )

# raw_effect > 0:
#   randomised matrices have higher Re(lambda_max) than the original web.
#
# standardised_effect < 0:
#   original web has lower Re(lambda_max) than the randomised null expectation.

#------------------------------------------------------------
# 11. Summary statistics
#------------------------------------------------------------

original_stats <- original_results %>%
  summarise(
    n_web = n(),
    median_Re_original = median(Re_original, na.rm = TRUE),
    q25_Re_original = quantile(Re_original, 0.25, na.rm = TRUE),
    q75_Re_original = quantile(Re_original, 0.75, na.rm = TRUE),
    mean_Re_original = mean(Re_original, na.rm = TRUE),
    sd_Re_original = sd(Re_original, na.rm = TRUE)
  )

rand_stats <- scenario_summary_by_web %>%
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

ses_stats <- effect_size_results %>%
  group_by(scenario) %>%
  summarise(
    n_web = n(),
    median_SES = median(standardised_effect, na.rm = TRUE),
    q25_SES = quantile(standardised_effect, 0.25, na.rm = TRUE),
    q75_SES = quantile(standardised_effect, 0.75, na.rm = TRUE),
    mean_SES = mean(standardised_effect, na.rm = TRUE),
    sd_SES = sd(standardised_effect, na.rm = TRUE),
    prop_negative = mean(standardised_effect < 0, na.rm = TRUE),
    prop_below_minus2 = mean(standardised_effect < -2, na.rm = TRUE),
    .groups = "drop"
  )

effect_tests <- effect_size_results %>%
  group_by(scenario) %>%
  summarise(
    n_web = n(),
    mean_raw_effect = mean(raw_effect, na.rm = TRUE),
    median_raw_effect = median(raw_effect, na.rm = TRUE),
    sd_raw_effect = sd(raw_effect, na.rm = TRUE),
    t_test_p = t.test(raw_effect, mu = 0)$p.value,
    wilcox_p = wilcox.test(raw_effect, mu = 0)$p.value,
    mean_standardised_effect = mean(standardised_effect, na.rm = TRUE),
    median_standardised_effect = median(standardised_effect, na.rm = TRUE),
    .groups = "drop"
  )

original_stats
rand_stats
ses_stats
effect_tests

#------------------------------------------------------------
# 12. Scenario comparison
#------------------------------------------------------------

scenario_means <- randomisation_results %>%
  group_by(
    valid_id,
    S_label,
    S_target,
    richness_eq,
    connectance_eq,
    maxTL_eq,
    scenario
  ) %>%
  summarise(
    mean_Re = mean(Re_lmax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    log_richness = log(richness_eq),
    log_connectance = log(connectance_eq)
  )

# Simple mixed model
lmm_fit <- lmer(
  mean_Re ~ scenario + log_richness + log_connectance + maxTL_eq + (1 | valid_id),
  data = scenario_means
)

summary(lmm_fit)
anova(lmm_fit)

emm <- emmeans(lmm_fit, ~ scenario)

pairwise_scenario_comparisons <- pairs(
  emm,
  adjust = "tukey"
)

pairwise_df <- as.data.frame(pairwise_scenario_comparisons)

pairwise_df

# Contrasts against no dependency
contrasts_vs_no <- contrast(
  emm,
  method = list(
    "full_dependency - no_dependency" = c(1, 0, 0, 0, -1),
    "c_dependency - no_dependency" = c(0, 1, 0, 0, -1),
    "B_dependency - no_dependency" = c(0, 0, 1, 0, -1),
    "e_dependency - no_dependency" = c(0, 0, 0, 1, -1)
  ),
  adjust = "holm"
)

contrasts_vs_no_df <- as.data.frame(contrasts_vs_no)

contrasts_vs_no_df

