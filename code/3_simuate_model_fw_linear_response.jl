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
S_bins = [
    (label = "S6_10", min = 6, max = 10, target = 8),
    (label = "S11_15", min = 11, max = 15, target = 13),
    (label = "S16_20", min = 16, max = 20, target = 18),
    (label = "S21_25", min = 21, max = 25, target = 23),
    (label = "S26_30", min = 26, max = 30, target = 28)
]

n_rep_per_Sbin = 100
n_valid_target = length(S_bins) * n_rep_per_Sbin

S_range = 10:100
C_range = (0.03, 0.40)
C_eq_range = (0.03, 0.35)

max_attempts = 200_000

alpha_range = (0.1, 5.0)
Z = 100
tmax = 10000
tmin = 1000

extinction_threshold = 1e-12
Re_ode_tolerance = 1e-6

out_dir = "data/model_fws/linear_response_sim"
mkpath(out_dir)

out_file = joinpath(out_dir, "model_food_webs.jld2")
feather_file = joinpath(out_dir, "model_food_webs.feather")

# ============================================================
# Functions
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
    open_bins = [
        S_bin for S_bin in S_bins
        if counts[S_bin.label] < n_rep
    ]

    isempty(open_bins) && return nothing

    S_bin = rand(open_bins)

    S_low = max(first(S_range), floor(Int, S_bin.target / 0.95))
    S_high = min(last(S_range), ceil(Int, S_bin.target / 0.55))

    S_low > S_high && return nothing

    S = rand(S_low:S_high)

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

        fw = Foodweb(
            :niche;
            S = S,
            C = C,
            reject_cycles = true,
            tol_C = 0.1
        )

        A = get_A(fw)
        sum(A .> 0) == 0 && continue

        # ------------------------------------------------------------
        # First simulation: burn-up on the initial generated web
        # ------------------------------------------------------------
        alpha, omega = random_linear_response_params(fw; alpha_range=alpha_range)

        model_burn = default_model(
            fw,
            LinearResponse(; alpha=alpha, w=omega),
            BodyMass(; Z=Z),
            Efficiency(:Miele2019)
        )

        B0 = rand(S) .+ 0.1

        sol_burn = simulate(
            model_burn,
            B0,
            tmax;
            callback=CallbackSet(
                extinction_callback(model_burn, extinction_threshold; verbose=false),
                TerminateSteadyState(
                    1e-8,
                    1e-6,
                    DiffEqCallbacks.allDerivPass;
                    min_t=tmin
                )
            ),
            show_degenerated=false
        )

        B_burn = Array(sol_burn.u[end])
        all(isfinite.(B_burn)) || continue

        # ------------------------------------------------------------
        # Prune burn-up result to realised alive-connected web
        # ------------------------------------------------------------
        A_burn0 = get_A(model_burn)
        is_producer_burn0 = vec(sum(A_burn0 .> 0, dims=2) .== 0)

        alive_burn = findall(B_burn .> extinction_threshold)
        isempty(alive_burn) && continue

        keep_burn = prune_realised_species(
            A_burn0,
            alive_burn,
            is_producer_burn0
        )

        isempty(keep_burn) && continue

        A_burn = A_burn0[keep_burn, keep_burn]
        sum(A_burn .> 0) == 0 && continue

        B0_finalrun = B_burn[keep_burn]

        alpha_burn = alpha[keep_burn]
        omega_burn = Matrix(omega)[keep_burn, keep_burn]

        # Rebuild the food web from the burn-up realised adjacency matrix
        fw_burn = Foodweb(A_burn)

        # ------------------------------------------------------------
        # Second simulation: final equilibrium on realised burn-up web
        # ------------------------------------------------------------
        model = default_model(
            fw_burn,
            LinearResponse(; alpha=alpha_burn, w=omega_burn),
            BodyMass(; Z=Z),
            Efficiency(:Miele2019)
        )

        sol = simulate(
            model,
            B0_finalrun,
            tmax;
            callback=CallbackSet(
                extinction_callback(model, extinction_threshold; verbose=false),
                TerminateSteadyState(
                    1e-8,
                    1e-6,
                    DiffEqCallbacks.allDerivPass;
                    min_t=tmin
                )
            ),
            show_degenerated=false
        )

        B_final = Array(sol.u[end])
       
        all(isfinite.(B_final)) || continue
        all(B_final .> extinction_threshold) || continue

        # Parameters from the rebuilt final model
        e, r, K, x = model.e, model.r, model.K, model.x

        eq = equilibrium_foodweb_summary(
            fw_burn,
            B_final,
            alpha_burn,
            omega_burn,
            e,
            r,
            K,
            x;
            extinction_threshold=extinction_threshold
        )

        isnothing(eq) && continue

        if eq.Re_ode_eq > Re_ode_tolerance
            fail_counts["unstable_ode"] = get(fail_counts, "unstable_ode", 0) + 1

            if fail_counts["unstable_ode"] <= 10
                println(
                    "Rejected unstable ODE case | ",
                    "S_eq = ", eq.S_eq,
                    " | C_eq = ", round(eq.C_eq, digits=3),
                    " | Re_ode = ", eq.Re_ode_eq
                )
            end

            continue
        end

        S_bin = match_S_bin(
            eq.S_eq,
            S_bins,
            target_counts,
            n_rep_per_Sbin
        )

        isnothing(S_bin) && continue

        if !(C_eq_range[1] <= eq.C_eq <= C_eq_range[2])
            continue
        end

        final_persistence = eq.S_eq / S

        global valid_id += 1
        target_counts[S_bin.label] += 1

        push!(
            valid_food_webs,
            (
                valid_id = valid_id,

                S_label = S_bin.label,
                S_min = S_bin.min,
                S_max = S_bin.max,
                S_target = S_bin.target,
                rep = target_counts[S_bin.label],

                initial_S = S,
                initial_C = C,
                initial_A = A,
                initial_L = sum(A .> 0),

                alpha = alpha,
                omega = omega,
                α = alpha,
                ω = omega,
                e = e,
                r = r,
                K = K,
                x = x,

                B_final = B_final,
                alive_persistence = persistence(B_final; threshold = extinction_threshold),
                final_persistence = final_persistence,

                keep_species = eq.keep_species,
                S_eq = eq.S_eq,
                L_eq = eq.L_eq,
                C_eq = eq.C_eq,
                maxTL_eq = eq.maxTL_eq,
                A_eq = eq.A_eq,
                B_eq = eq.B_eq,
                TL_eq = eq.TL_eq,
                is_producer_eq = eq.is_producer_eq,

                alpha_eq = eq.alpha_eq,
                omega_eq = eq.omega_eq,
                α_eq = eq.α_eq,
                ω_eq = eq.ω_eq,
                e_eq = eq.e_eq,
                r_eq = eq.r_eq,
                K_eq = eq.K_eq,
                x_eq = eq.x_eq,

                J_emp = eq.J_emp_eq,
                J_ode = eq.J_ode_eq,
                Re_emp = eq.Re_emp_eq,
                Re_ode = eq.Re_ode_eq,
                positive_diag_idx = eq.positive_diag_idx,
                positive_diag_values = eq.positive_diag_values,

                totalB_final = sum(B_final),
                totalB_eq = sum(eq.B_eq)
            )
        )

        println(
            "Accepted $valid_id / $n_valid_target | ",
            "S_bin = ", S_bin.label,
            " | rep = ", target_counts[S_bin.label],
            " | S_eq = ", eq.S_eq,
            " | C_eq = ", round(eq.C_eq, digits = 3),
            " | maxTL_eq = ", round(eq.maxTL_eq, digits = 3),
            " | Re_emp = ", round(eq.Re_emp_eq, digits = 3),
            " | Re_ode = ", round(eq.Re_ode_eq, digits = 4)
        )

        if attempt % 100 == 0
            n_accepted = sum(values(target_counts))
            n_failed   = sum(values(fail_counts))
            n_tracked  = n_accepted + n_failed
            n_missing  = attempt - n_tracked
            
            println("Attempt $attempt | accepted $n_accepted / $n_valid_target")
            println("Counts: ", target_counts)
            println("Failures: ", fail_counts)
            println("Tracked: $n_tracked / $attempt | Missing: $n_missing")
        end
    
    catch err
        msg = sprint(showerror, err)

        fail_type =
            occursin("Connectance", msg) ? "connectance" :
            occursin("NaN", msg) || occursin("Inf", msg) ? "nan_or_inf" :
            "other"

        fail_counts[fail_type] = get(fail_counts, fail_type, 0) + 1

        if attempt % 100 == 0
            n_accepted = sum(values(target_counts))
            n_failed   = sum(values(fail_counts))
            n_tracked  = n_accepted + n_failed
            n_missing  = attempt - n_tracked
            
            println("Attempt $attempt | accepted $n_accepted / $n_valid_target")
            println("Counts: ", target_counts)
            println("Failures: ", fail_counts)
            println("Tracked: $n_tracked / $attempt | Missing: $n_missing")
        end
    end
end


#=
# check
diag_df = DataFrame(map(valid_food_webs) do fw
    J = fw.J_ode

    diag_vals = abs.(diag(J))

    Joff = copy(J)
    Joff[diagind(Joff)] .= 0.0
    off_vals = abs.(Joff[Joff .!= 0])

    (
        valid_id = fw.valid_id,
        S_eq = fw.S_eq,
        C_eq = fw.C_eq,
        mean_diag_over_offdiag = mean(diag_vals) / mean(off_vals),
        median_diag_over_offdiag = median(diag_vals) / median(off_vals),
        max_diag_over_offdiag = maximum(diag_vals) / maximum(off_vals)
    )
end)
describe(diag_df)
=#

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
    C_eq_range = C_eq_range,
    Re_ode_tolerance = Re_ode_tolerance
)

# ============================================================
# Save R Feather file
# ============================================================
dep_rows = DataFrame[]

for fw in valid_food_webs
    richness_eq = fw.S_eq
    connectance_eq = fw.C_eq
    maxTL_eq = fw.maxTL_eq

    alpha_mat = reshape(fw.alpha_eq, :, 1)
    cij_eq = alpha_mat .* Matrix(fw.omega_eq)

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
    J_ode_diag_prey = Float64[]
    J_ode_diag_predator = Float64[]

    for consumer in 1:richness_eq, resource in 1:richness_eq
        if cij_eq[consumer, resource] > 0
            push!(prey_idx, resource)
            push!(predator_idx, consumer)

            push!(cij_value, cij_eq[consumer, resource])
            push!(alpha_value, fw.alpha_eq[consumer])
            push!(omega_value, fw.omega_eq[consumer, resource])

            push!(B_prey, fw.B_eq[resource])
            push!(B_predator, fw.B_eq[consumer])
            push!(e_value, fw.e_eq[consumer, resource])

            push!(J_topdown, fw.J_emp[resource, consumer])
            push!(J_bottomup, fw.J_emp[consumer, resource])
            push!(J_ode_diag_prey, fw.J_ode[resource, resource])
            push!(J_ode_diag_predator, fw.J_ode[consumer, consumer])
        end
    end

    n_links = length(cij_value)
    n_links == 0 && continue

    push!(
        dep_rows,
        DataFrame(
            valid_id = fill(fw.valid_id, n_links),

            S_label = fill(fw.S_label, n_links),
            S_min = fill(fw.S_min, n_links),
            S_max = fill(fw.S_max, n_links),
            S_target = fill(fw.S_target, n_links),
            rep = fill(fw.rep, n_links),

            richness_eq = fill(richness_eq, n_links),
            connectance_eq = fill(connectance_eq, n_links),
            maxTL_eq = fill(maxTL_eq, n_links),
            L_eq = fill(fw.L_eq, n_links),

            Re_emp = fill(fw.Re_emp, n_links),
            Re_ode = fill(fw.Re_ode, n_links),

            prey = prey_idx,
            predator = predator_idx,
            cij = cij_value,
            alpha = alpha_value,
            omega = omega_value,
            B_prey = B_prey,
            B_predator = B_predator,
            e = e_value,
            J_topdown = J_topdown,
            J_bottomup = J_bottomup,
            J_ode_diag_prey = J_ode_diag_prey,
            J_ode_diag_predator = J_ode_diag_predator
        )
    )
end

model_foodwebs_df = vcat(dep_rows...)
Arrow.write(feather_file, model_foodwebs_df)

