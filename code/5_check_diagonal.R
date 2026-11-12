# Diagonal and off-diagonal diagnostics for model food webs

rm(list = ls())

library(tidyverse)
library(here)
library(arrow)

source(here("code/src/internal.R"))

feather_file <- here("data/model_fws/model_food_webs_dependency_links.feather")
reference_suffix <- "_Re_ode"
S_label_levels <- c("6_10", "11_15", "16_20", "21_25", "26_30", "31_35", "36_40")

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

web_meta <- modelfw_df %>%
  distinct(valid_id, S_label, richness_eq, connectance_eq, maxTL_eq, Re_ode)

dir.create(here("figures"), recursive = TRUE, showWarnings = FALSE)

theme_dep <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
      axis.line = element_blank(),
      axis.text = element_text(colour = "black", size = base_size - 2),
      axis.title = element_text(colour = "black", size = base_size),
      axis.ticks = element_line(linewidth = 0.45),
      legend.position = "bottom",
      legend.text = element_text(size = base_size - 4),
      plot.tag = element_text(face = "bold", size = base_size)
    )
}

jacobians <- map(web_objects, function(web) {
  build_jacobian(web$C, web$B, web$e_link, J_diag = web$J_ode_diag)
})

# C[i, j] is the interaction coefficient for consumer j feeding on resource i.
jacobian_diagonal_elements <- map_dfr(names(web_objects), function(id) {
  
  C <- web_objects[[id]]$C
  diagonal <- diag(jacobians[[id]])
  row_links <- rowSums(C != 0)  # number of predators
  col_links <- colSums(C != 0)  # number of resources
  
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
  left_join(web_meta, by = "valid_id")

Jii_role_summary <- jacobian_diagonal_elements %>%
  group_by(trophic_role) %>%
  summarise(
    n = n(),
    median_abs_Jii = median(abs_Jii, na.rm = TRUE),
    q25 = quantile(abs_Jii, 0.25, na.rm = TRUE),
    q75 = quantile(abs_Jii, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(Jii_role_summary)

web_trophic_composition <- jacobian_diagonal_elements %>%
  group_by(valid_id) %>%
  summarise(
    n_species = n(),
    n_basal = sum(trophic_role == "Basal"),
    n_intermediate = sum(trophic_role == "Intermediate consumer"),
    n_top = sum(trophic_role == "Top consumer"),
    prop_basal = mean(trophic_role == "Basal"),
    prop_intermediate = mean(trophic_role == "Intermediate consumer"),
    prop_top = mean(trophic_role == "Top consumer"),
    .groups = "drop"
  )

jacobian_diagnostics <- map_dfr(jacobians, function(J) {
  off_diagonal <- J[row(J) != col(J) & J != 0]
  
  tibble(
    median_abs_diagonal = median(abs(diag(J)), na.rm = TRUE),
    median_abs_off_diagonal = median(abs(off_diagonal), na.rm = TRUE),
    median_abs_positive_off_diagonal =
      median(off_diagonal[off_diagonal > 0], na.rm = TRUE),
    median_abs_negative_off_diagonal =
      median(-off_diagonal[off_diagonal < 0], na.rm = TRUE)
  )
}, .id = "valid_id") %>%
  left_join(web_meta, by = "valid_id") %>%
  left_join(web_trophic_composition, by = "valid_id")

print(summary(jacobian_diagnostics))

# Plots
p_Jii_role <- jacobian_diagonal_elements %>%
  filter(is.finite(abs_Jii), abs_Jii > 0, trophic_role != "Isolated") %>%
  ggplot(aes(abs_Jii, fill = trophic_role)) +
  geom_histogram(bins = 80, position = "identity", alpha = 0.55) +
  scale_x_log10() +
  labs(
    x = expression(abs(J[ii])),
    y = "Number of species",
    fill = "Trophic role"
  ) +
  theme_dep(20)

ggsave(
  here("figures", "FigureS2_left.jpeg"),
  p_Jii_role, width = 10, height = 7, dpi = 600
)


p_diag_offdiag <- ggplot(
  jacobian_diagnostics,
  aes(median_abs_diagonal, median_abs_off_diagonal, colour = prop_basal)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_viridis_c(
    option = "cividis",
    limits = range(jacobian_diagnostics$prop_basal, na.rm = TRUE),
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

ggsave(
  here("figures", "FigureS2_right.jpeg"),
  p_diag_offdiag, width = 8, height = 7, dpi = 600
)
