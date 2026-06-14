
using Pkg
Pkg.activate(".")

using EcologicalNetworksDynamics
using Random
using Distributions
using LinearAlgebra
using SciMLBase: CallbackSet
using DiffEqCallbacks
using SparseArrays
using JLD2
using DataFrames
using Arrow

include("internal.jl")

Random.seed!(2026)

# ============================================================
# Settings
# ============================================================
# Realised richness bins.
# Connectance is not binned; it is kept as a continuous realised variable.
S_bins = [
    (label = "S6_10", min = 6, max = 10, target = 8),
    (label = "S11_15", min = 11, max = 15, target = 13),
    (label = "S16_20", min = 16, max = 20, target = 18),
    (label = "S21_25", min = 21, max = 25, target = 23),
    (label = "S26_30", min = 26, max = 30, target = 28),
    (label = "S31_35", min = 31, max = 35, target = 33),
    (label = "S36_40", min = 36, max = 40, target = 38)
]

n_rep_per_Sbin = 5
n_valid_target = length(S_bins) * n_rep_per_Sbin

# Initial niche-model ranges.
# These are candidate-generation ranges only.
S_range = 10:90
C_range = (0.03, 0.40)

# Keep realised connectance within a broad range.
# This avoids extremely sparse or extremely dense realised webs,
# but does not force connectance into bins.
C_eq_range = (0.03, 0.35)

max_attempts = 200_000

α_range = (0.1, 5.0)
Z = 100
tmax = 2000
tmin = 500

extinction_threshold = 1e-12

out_dir = "data/model_fws/linear_response_niche_webs"
mkpath(out_dir)

out_file = joinpath(out_dir, "valid_food_webs_varying_SC.jld2")
feather_file = joinpath(out_dir, "valid_food_webs_varying_SC.feather")

# ============================================================
# Helper functions
# ============================================================
function empty_S_counts(S_bins)
    Dict(S_bin.label => 0 for S_bin in S_bins)
end

function match_S_bin(S_eq, S_bins, counts, n_rep)
    for S_bin in S_bins
        in_S = S_bin.min <= S_eq <= S_bin.max

        if in_S && counts[S_bin.label] < n_rep
            return S_bin
        end
    end

    return nothing
end

function propose_initial_web(S_bins, counts, n_rep, S_range, C_range)
    # Choose an unfilled realised-richness bin.
    open_bins = [
        S_bin for S_bin in S_bins
        if counts[S_bin.label] < n_rep
    ]

    isempty(open_bins) && return nothing

    S_bin = rand(open_bins)

    # Choose initial S somewhat larger than the realised target,
    # because extinctions and disconnected-species removal may reduce richness.
    S_low = max(first(S_range), floor(Int, S_bin.target / 0.95))
    S_high = min(last(S_range), ceil(Int, S_bin.target / 0.55))

    S_low > S_high && return nothing

    S = rand(S_low:S_high)

    # Choose initial C broadly.
    # Also respect the minimum connectance allowed by the niche generator.
    C_min_allowed = (S - 1) / S^2
    C_low = max(C_range[1], C_min_allowed)
    C_high = C_range[2]

    C_low >= C_high && return nothing

    C = rand(Uniform(C_low, C_high))

    return S, C, S_bin.label
end

# ============================================================
# Run
# ============================================================
target_counts = empty_S_counts(S_bins)
fail_counts = Dict{String, Int}()

valid_food_webs = []

attempt = 0
valid_id = 0

while valid_id < n_valid_target && attempt < max_attempts

    global attempt += 1

    try
        proposed = propose_initial_web(
            S_bins,
            target_counts,
            n_rep_per_Sbin,
            S_range,
            C_range
        )

        isnothing(proposed) && continue

        S, C, proposed_S_bin = proposed

        # ----- 1. Generate niche model food web -----

        fw = Foodweb(
            :niche;
            S = S,
            C = C,
            reject_cycles = true,
            tol_C = 0.1
        )

        A = get_A(fw)
        sum(A .> 0) == 0 && continue


        # ----- 2. Assign parameter values -----
        α, ω = random_linear_response_params(fw; α_range = α_range)

        model = default_model(
            fw,
            LinearResponse(; alpha = α, w = ω),
            BodyMass(; Z = Z),
            Efficiency(:Miele2019)
        )

        B0 = rand(S) .+ 0.1

        # ----- 3. Simulate dynamics -----
        sol = simulate(
            model,
            B0,
            tmax;
            callback = CallbackSet(
                extinction_callback(model, extinction_threshold; verbose = false),
                TerminateSteadyState(
                    1e-8,
                    1e-6,
                    DiffEqCallbacks.allDerivPass;
                    min_t = tmin
                )
            ),
            show_degenerated = false
        )

        B_final = Array(sol.u[end])
        all(isfinite.(B_final)) || continue

        # ----- 4. Equilibrium food web, Jacobian matrix and stability -----
        e, r, K, x = model.e, model.r, model.K, model.x

        eq = equilibrium_foodweb_summary(
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

        isnothing(eq) && continue

        # Keep realised richness from 10 to 50.
        S_bin = match_S_bin(
            eq.S_eq,
            S_bins,
            target_counts,
            n_rep_per_Sbin
        )

        isnothing(S_bin) && continue

        # Keep realised connectance varying, but within a broad valid range.
        if !(C_eq_range[1] <= eq.C_eq <= C_eq_range[2])
            continue
        end

        final_persistence = eq.S_eq / S

        # ----- 5. Save metadata and results -----
        global valid_id += 1
        target_counts[S_bin.label] += 1

        push!(
            valid_food_webs,
            (
                valid_id = valid_id,

                # realised richness-bin information
                S_label = S_bin.label,
                S_min = S_bin.min,
                S_max = S_bin.max,
                S_target = S_bin.target,
                rep = target_counts[S_bin.label],

                # initial web
                initial_S = S,
                initial_C = C,
                initial_A = A,
                initial_L = sum(A .> 0),

                # original parameters
                α = α,
                ω = ω,
                e = e,
                r = r,
                K = K,
                x = x,

                # dynamics
                B_final = B_final,
                alive_persistence = persistence(B_final; threshold = extinction_threshold),
                final_persistence = final_persistence,

                # realised equilibrium web
                keep_species = eq.keep_species,
                S_eq = eq.S_eq,
                L_eq = eq.L_eq,
                C_eq = eq.C_eq,
                maxTL_eq = eq.maxTL_eq,
                A_eq = eq.A_eq,
                B_eq = eq.B_eq,
                TL_eq = eq.TL_eq,

                # realised parameters
                α_eq = eq.α_eq,
                ω_eq = eq.ω_eq,
                e_eq = eq.e_eq,
                r_eq = eq.r_eq,
                K_eq = eq.K_eq,
                x_eq = eq.x_eq,

                # Jacobians and stability
                J_emp = eq.J_emp_eq,
                J_ode = eq.J_ode_eq,
                Re_emp = eq.Re_emp_eq,
                Re_ode = eq.Re_ode_eq,

                # biomass summaries
                totalB_final = sum(B_final),
                totalB_eq = sum(eq.B_eq)
            )
        )

        println(
            "Accepted $valid_id / $n_valid_target | ",
            "S_bin = ", S_bin.label,
            " | rep = ", target_counts[S_bin.label],
           # " | initial S = ", S,
           # " | initial C = ", round(C, digits = 3),
            " | S_eq = ", eq.S_eq,
            " | C_eq = ", round(eq.C_eq, digits = 3),
            " | maxTL_eq = ", round(eq.maxTL_eq, digits = 3),
            " | Re_emp = ", round(eq.Re_emp_eq, digits = 3),
            " | Re_ode = ", round(eq.Re_ode_eq, digits = 3)
        )

        if attempt % 100 == 0
            println("Attempt $attempt | accepted $valid_id / $n_valid_target")
            println("Counts: ", target_counts)
            println("Failures: ", fail_counts)
        end

    catch err
        msg = sprint(showerror, err)

        fail_type =
            occursin("Connectance", msg) ? "connectance" :
            occursin("NaN", msg) || occursin("Inf", msg) ? "nan_or_inf" :
            "other"

        fail_counts[fail_type] = get(fail_counts, fail_type, 0) + 1

        if attempt % 100 == 0
            println("Attempt $attempt | accepted $valid_id / $n_valid_target")
            println("Counts: ", target_counts)
            println("Failures: ", fail_counts)
        end
    end
end

# quick check
valid_food_webs[1].Re_emp

# ============================================================
# Save full Julia output
# ============================================================
jldsave(
    out_file;
    valid_food_webs = valid_food_webs,
    S_bins = S_bins,
    target_counts = target_counts,
    fail_counts = fail_counts,
    n_rep_per_Sbin = n_rep_per_Sbin,
    C_eq_range = C_eq_range
)

# ============================================================
# Save R Feather file
# ============================================================
dep_rows = DataFrame[]

for fw in valid_food_webs
    richness_eq = fw.S_eq
    connectance_eq = fw.C_eq
    maxTL_eq = fw.maxTL_eq

    # c_{consumer, resource} = α_consumer * ω_{consumer, resource}
    α_mat = reshape(fw.α_eq, :, 1)
    cij_eq = α_mat .* Matrix(fw.ω_eq)

    prey_idx = Int[]
    predator_idx = Int[]
    cij_value = Float64[]
    alpha_value = Float64[]
    omega_value = Float64[]
    B_prey = Float64[]
    B_predator = Float64[]
    e_value = Float64[]
    J_topdown = Float64[]
    J_bottomup = Float64[]

    for consumer in 1:richness_eq, resource in 1:richness_eq
        if cij_eq[consumer, resource] > 0
            push!(prey_idx, resource)
            push!(predator_idx, consumer)

            push!(cij_value, cij_eq[consumer, resource])
            push!(alpha_value, fw.α_eq[consumer])
            push!(omega_value, fw.ω_eq[consumer, resource])

            push!(B_prey, fw.B_eq[resource])
            push!(B_predator, fw.B_eq[consumer])
            push!(e_value, fw.e_eq[consumer, resource])

            # J_emp[resource, consumer]: top-down effect of predator on prey
            # J_emp[consumer, resource]: bottom-up effect of prey on predator
            push!(J_topdown, fw.J_emp[resource, consumer])
            push!(J_bottomup, fw.J_emp[consumer, resource])
        end
    end

    n_links = length(cij_value)

    if n_links == 0
        continue
    end

    push!(
        dep_rows,
        DataFrame(
            valid_id = fill(fw.valid_id, n_links),

            # realised richness-bin information
            S_label = fill(fw.S_label, n_links),
            S_min = fill(fw.S_min, n_links),
            S_max = fill(fw.S_max, n_links),
            S_target = fill(fw.S_target, n_links),
            rep = fill(fw.rep, n_links),

            # realised food-web summary
            richness_eq = fill(richness_eq, n_links),
            connectance_eq = fill(connectance_eq, n_links),
            maxTL_eq = fill(maxTL_eq, n_links),
            L_eq = fill(fw.L_eq, n_links),

            # stability
            Re_emp = fill(fw.Re_emp, n_links),
            Re_ode = fill(fw.Re_ode, n_links),

            # link-level data
            prey = prey_idx,
            predator = predator_idx,
            cij = cij_value,
            alpha = alpha_value,
            omega = omega_value,
            B_prey = B_prey,
            B_predator = B_predator,
            e = e_value,
            J_topdown = J_topdown,
            J_bottomup = J_bottomup
        )
    )
end

niche_foodwebs_df = vcat(dep_rows...)
Arrow.write(feather_file, niche_foodwebs_df)

