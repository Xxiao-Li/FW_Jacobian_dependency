# 

This repository contains the code and data used in the paper "Structural dependencies in trophic interaction strengths explain stability in natural food webs".

In this study, we tested how shared biological quantities link trophic interaction strengths and affect local stability in a three-species intraguild predation module, Lake Constance food webs and bioenergetic model food webs. Retaining all structural dependencies reduced the destabilising effects of interspecific interactions most strongly across all three systems, with larger effects in richer and more connected model food webs.

# Analysis workflow

Open `Jacobian_dependency.Rproj` in RStudio or run the scripts from the repository root. The scripts read from `data/` and save plots to `figures/`.

-   [\`1_three_species_modules.R\`](code/1_three_species_modules.R) compares dependency scenarios in the three-species module (Figure 2).

-   [\`2_LC_food_webs.R\`](code/2_LC_food_webs.R) compares empirical and randomised Lake Constance food webs (Figure 3).

-   [\`4_analyse_model_fws.R\`](code/4_analyse_model_fws.R) analyses the model food webs using the saved simulation and randomisation results (Figure 4 by default). For Figure 5, set `reference_type <- "Re_ode"` to use the full Jacobian and switch the output filenames to `Figure5`.

-   [\`5_check_diagonal.R\`](code/5_check_diagonal.R) checks the diagonal and off-diagonal entries of the model Jacobians (Figure S1).

### Generate model food webs in Julia

1.  Run [\`3_1_HPC_simuate_model_fws.jl\`](code/3_1_HPC_simuate_model_fws.jl) for array tasks 1–70 to generate the model food webs.

2.  Run [\`3_2_combine_array_outputs.jl\`](code/3_2_combine_array_outputs.jl) to combine the task outputs into one Feather file, then place it in `data/model_fws/`.
