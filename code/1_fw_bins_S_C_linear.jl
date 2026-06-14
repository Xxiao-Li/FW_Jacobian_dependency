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

include("internal.jl")

Random.seed!(2026)

# ============================================================
# Settings
# ============================================================
S_bins = [
    (label = "low",    target = 10, min = 8,  max = 12),
    (label = "medium", target = 20, min = 16, max = 24),
    (label = "high",   target = 40, min = 32, max = 48)
]

C_bins = [
    (label = "low",    target = 0.10, min = 0.09, max = 0.12),
    (label = "medium", target = 0.17, min = 0.15, max = 0.19),
    (label = "high",   target = 0.26, min = 0.23, max = 0.29)
]

n_rep_per_combo = 1
n_valid_target = length(S_bins) * length(C_bins) * n_rep_per_combo

S_range = 10:80
C_range = (0.03, 0.40)

max_attempts = 100_000

α_range = (1.0, 10.0)
Z = 100
tmax = 3_000
tmin = 1_000

extinction_threshold = 1e-8

out_dir = "data/model_fws/linear_response_niche_webs"
mkpath(out_dir)

out_file = joinpath(out_dir, "valid_food_webs_realised_SC_bins.jld2")

# ============================================================
# Helper functions
# ============================================================
function empty_counts(S_bins, C_bins)
    Dict((S.label, C.label) => 0 for S in S_bins, C in C_bins)
end

function match_bin(S_eq, C_eq, S_bins, C_bins, counts, n_rep)
    for S_bin in S_bins, C_bin in C_bins
        in_S = S_bin.min <= S_eq <= S_bin.max
        in_C = C_bin.min <= C_eq <= C_bin.max
        key = (S_bin.label, C_bin.label)

        if in_S && in_C && counts[key] < n_rep
            return S_bin, C_bin, key
        end
    end

    return nothing
end

function propose_initial_web(S_bins, C_bins, counts, n_rep, S_range, C_range)
    # Choose an unfilled bin
    open_bins = [
        (S_bin, C_bin) for S_bin in S_bins, C_bin in C_bins
        if counts[(S_bin.label, C_bin.label)] < n_rep
    ]

    isempty(open_bins) && return nothing

    S_bin, C_bin = rand(open_bins)

    # Choose initial S slightly larger than target because extinctions may occur
    S_low = max(first(S_range), floor(Int, S_bin.target / 0.95))
    S_high = min(last(S_range), ceil(Int, S_bin.target / 0.55))
    S = rand(S_low:S_high)

    # Respect EcologicalNetworksDynamics minimum connectance
    C_min_allowed = (S - 1) / S^2
    C_low = max(C_range[1], C_min_allowed, C_bin.target - 0.06)
    C_high = min(C_range[2], C_bin.target + 0.08)

    C_low >= C_high && return nothing

    C = rand(Uniform(C_low, C_high))

    return S, C
end

# ============================================================
# Run
# ============================================================
target_counts = empty_counts(S_bins, C_bins)
fail_counts = Dict{String, Int}()

valid_food_webs = []

attempt = 0
valid_id = 0

while valid_id < n_valid_target && attempt < max_attempts

    global attempt += 1

    try
        proposed = propose_initial_web(
            S_bins,
            C_bins,
            target_counts,
            n_rep_per_combo,
            S_range,
            C_range
        )

        isnothing(proposed) && continue

        S, C = proposed

        fw = Foodweb(
            :niche;
            S = S,
            C = C,
            reject_cycles = true,
            tol_C = 0.1
        )

        A = get_A(fw)
        sum(A .> 0) == 0 && continue

        α, ω = random_linear_response_params(fw; α_range = α_range)

        model = default_model(
            fw,
            LinearResponse(; alpha = α, w = ω),
            BodyMass(; Z = Z),
            Efficiency(:Miele2019)
        )

        B0 = rand(S) .+ 0.1

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

        final_persistence = eq.S_eq / S

        matched = match_bin(
            eq.S_eq,
            eq.C_eq,
            S_bins,
            C_bins,
            target_counts,
            n_rep_per_combo
        )

        isnothing(matched) && continue

        S_bin, C_bin, key = matched

        global valid_id += 1
        target_counts[key] += 1

        push!(
            valid_food_webs,
            (
                valid_id = valid_id,

                # bin information
                S_label = S_bin.label,
                C_label = C_bin.label,
                S_target = S_bin.target,
                C_target = C_bin.target,
                rep = target_counts[key],

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
                A_eq = eq.A_eq,
                B_eq = eq.B_eq,

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
            " | C_bin = ", C_bin.label,
            " | rep = ", target_counts[key],
            " | initial S = ", S,
            " | initial C = ", round(C, digits = 3),
            " | S_eq = ", eq.S_eq,
            " | C_eq = ", round(eq.C_eq, digits = 3)
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


# ============================================================
# Save
# ============================================================
jldsave(
    out_file;
    valid_food_webs = valid_food_webs,
    S_bins = S_bins,
    C_bins = C_bins,
    target_counts = target_counts,
    fail_counts = fail_counts,
    n_rep_per_combo = n_rep_per_combo
)

println("Finished.")
println("Valid food webs: ", length(valid_food_webs))
println("Attempts: ", attempt)
println("Target counts: ", target_counts)
println("Failure counts: ", fail_counts)
println("Saved to: ", out_file)