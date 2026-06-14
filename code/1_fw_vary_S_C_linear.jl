using Pkg
Pkg.activate(".")

using EcologicalNetworksDynamics
using Random
using Distributions
using SparseArrays
using Plots
using LinearAlgebra
using DiffEqCallbacks
using JLD2

include("internal.jl")

Random.seed!(2026)

n_valid_target = 10
max_attempts = 10_000

S_range = [10, 20, 30]
C_range = (0.05, 0.30)

α_range = (1.0, 10.0)
Z = 100
tmax = 5_000
tmin = 1_000

extinction_threshold = 1e-8
persistence_threshold = 0.6

out_dir = "data/model_fws/linear_response_niche_webs"
mkpath(out_dir)

# ============================================================
# Run until 100 valid food webs are obtained
# ============================================================

valid_food_webs = []

attempt = 0
valid_id = 0

while valid_id < n_valid_target && attempt < max_attempts

    global attempt += 1

    try

        # ----- 1. Generate niche model food webs -----
        S = rand(S_range)
        C = rand(Uniform(C_range[1], C_range[2]))

        fw = Foodweb(:niche; S = S, C = C)
        A = get_A(fw)

        # Skip empty or invalid topologies
        if sum(A .> 0) == 0
            continue
        end

        # ----- 2. Assign parameter values -----
        α, ω = random_linear_response_params(
            fw;
            α_range = α_range
        )

        model = default_model(
            fw,
            LinearResponse(; alpha = α, w = ω),
            BodyMass(; Z = Z),
            Efficiency(:Miele2019)
        )

        B0 = rand(S) .+ 0.1

        # ----- 3. Simulate the dynamics -----
        sol = simulate(
            model,
            B0,
            tmax;
            callback = CallbackSet(
                extinction_callback(model, extinction_threshold; verbose = false),
                TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass; min_t = tmin)
            ),
            show_degenerated = false
        )

        B_final = Array(sol.u[end])

        if !all(isfinite.(B_final))
            continue
        end

        alive_persistence = persistence(
            B_final;
            threshold = extinction_threshold
        )

        # First filter: alive species only
        if alive_persistence <= persistence_threshold
            continue
        end

        # ----- 4. Equilibrium food web, Jacobian matrix and stability -----
        e = model.e
        r = model.r
        K = model.K
        x = model.x

        eq_summary = equilibrium_foodweb_summary(
            fw,
            B_final,
            α,
            ω,
            e,
            r,
            K,
            x;
            extinction_threshold = extinction_threshold
        )

        if isnothing(eq_summary)
            continue
        end

        # Second filter: alive and connected species
        final_persistence = eq_summary.S_eq / S

        if final_persistence <= persistence_threshold
            continue
        end

        # ----- 5. Save meta data and results -----
        global valid_id += 1

        push!(
            valid_food_webs,
            (
                valid_id = valid_id,
                attempt = attempt,

                # initial food web
                initial_S = S,
                initial_C = C,
                initial_A = A,

                # original parameters
                α = α,
                ω = ω,
                e = e,
                r = r,
                K = K,
                x = x,

                # equilibrium realised food web
                keep_species = eq_summary.keep_species,

                S_eq = eq_summary.S_eq,
                L_eq = eq_summary.L_eq,
                C_eq = eq_summary.C_eq,
                A_eq = eq_summary.A_eq,
                B_eq = eq_summary.B_eq,

                # equilibrium parameters
                α_eq = eq_summary.α_eq,
                ω_eq = eq_summary.ω_eq,
                e_eq = eq_summary.e_eq,
                r_eq = eq_summary.r_eq,
                K_eq = eq_summary.K_eq,
                x_eq = eq_summary.x_eq,

                # equilibrium Jacobians
                J_emp = eq_summary.J_emp_eq,
                J_ode = eq_summary.J_ode_eq,
                Re_emp = eq_summary.Re_emp_eq,
                Re_ode = eq_summary.Re_ode_eq,
            )
        )

        println(
            "Accepted $valid_id / $n_valid_target | ",
            "S_initial = $S | ",
            "S_eq = ", eq_summary.S_eq, " | ",
            "C = ", round(C, digits = 3), " | ",
            "C_eq = ", round(eq_summary.C_eq, digits = 3), " | "
        )

    catch err
        println("Attempt $attempt failed: ", err)
        continue
    end
end


# ============================================================
# Save output
# ============================================================
jldsave(
    joinpath(out_dir, "valid_food_webs.jld2");
    valid_food_webs = valid_food_webs
)