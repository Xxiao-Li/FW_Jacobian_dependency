using JLD2
using DataFrames
using Arrow

input_dir = "data/model_fws/linear_response_sim_7bins_70tasks"

output_file = joinpath(
    input_dir,
    "model_food_webs_dependency_links.feather"
)

# Final task files only; checkpoint files are excluded
files = filter(
    f -> occursin(r"^model_food_webs_task_\d+\.jld2$", basename(f)),
    readdir(input_dir; join=true)
)

sort!(files)

# Read all accepted food webs
valid_food_webs = Any[]

for file in files
    append!(
        valid_food_webs,
        load(file, "valid_food_webs")
    )
end

println("Loaded $(length(valid_food_webs)) food webs")

ids = [fw.valid_id for fw in valid_food_webs]
missing_ids = setdiff(1:700, ids)
println("Missing IDs: ", missing_ids)


rows = DataFrame[]

for fw in valid_food_webs
    S = fw.S_eq

    alpha = collect(fw.alpha_eq)
    omega = Matrix(fw.omega_eq)
    e = Matrix(fw.e_eq)

    A = Matrix(fw.A_eq)
    B = collect(fw.B_eq)

    J_emp = Matrix(fw.J_emp)
    J_ode = Matrix(fw.J_ode)

    cij = reshape(alpha, :, 1) .* omega

    for predator in 1:S
        for prey in 1:S
            A[predator, prey] == 0 && continue

            push!(
                rows,
                DataFrame(
                    valid_id=fw.valid_id,
                    task_id=fw.task_id,
                    S_label=fw.S_label,
                    rep=fw.rep,

                    richness_eq=fw.S_eq,
                    connectance_eq=fw.C_eq,
                    maxTL_eq=fw.maxTL_eq,
                    L_eq=fw.L_eq,

                    Re_emp=fw.Re_emp,
                    Re_ode=fw.Re_ode,

                    prey=prey,
                    predator=predator,

                    cij=cij[predator, prey],
                    alpha=alpha[predator],
                    omega=omega[predator, prey],
                    e=e[predator, prey],

                    B_prey=B[prey],
                    B_predator=B[predator],

                    J_topdown=J_emp[prey, predator],
                    J_bottomup=J_emp[predator, prey],

                    J_ode_topdown=J_ode[prey, predator],
                    J_ode_bottomup=J_ode[predator, prey],

                    J_ode_diag_prey=J_ode[prey, prey],
                    J_ode_diag_predator=J_ode[predator, predator]
                )
            )
        end
    end
end

model_foodwebs_df = vcat(rows...)

sort!(
    model_foodwebs_df,
    [:valid_id, :predator, :prey]
)

Arrow.write(
    output_file,
    model_foodwebs_df
)
