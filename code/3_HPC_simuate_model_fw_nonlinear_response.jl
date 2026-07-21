#=
classic_response_h2_array_7bins_70tasks.jl

Nonlinear functional response:
- Hill exponent h = 2
- Intraspecific consumer interference i0 = 0
- This corresponds to a Type III functional response
=#

using Pkg
Pkg.activate(".")

using EcologicalNetworksDynamics
using Random
using Distributions
using LinearAlgebra
using SparseArrays
using StatsBase
using SciMLBase: CallbackSet
using DiffEqCallbacks
using JLD2

include("src/internal.jl")

# ============================================================
# Array settings
# ============================================================

length(ARGS) >= 1 ||
    error("Run as: julia script.jl TASK_ID [N_TASKS]")

task_id = parse(Int, ARGS[1])
n_tasks = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 70

1 <= task_id <= n_tasks ||
    error("task_id must be between 1 and $n_tasks")

base_seed = 2026
task_tag = lpad(string(task_id), length(string(n_tasks)), '0')

# ============================================================
# Simulation settings
# ============================================================

hill_exponent = 2.0
i0_intra = 0.0

# Bins are based on realised richness after burn-in,
# pruning and the final equilibrium simulation.
S_bins = [
    (label="S6_10",  min=6,  max=10, target=8),
    (label="S11_15", min=11, max=15, target=13),
    (label="S16_20", min=16, max=20, target=18),
    (label="S21_25", min=21, max=25, target=23),
    (label="S26_30", min=26, max=30, target=28),
    (label="S31_35", min=31, max=35, target=33),
    (label="S36_40", min=36, max=40, target=38)
]

# Number of accepted food webs across all tasks for each bin
n_rep_per_Sbin = 100

# Maximum attempts during one execution of each array task
max_attempts_this_run = 500_000

# Possible initial richness
S_range = 10:150

# Initial connectance proposal range
C_range = (0.03, 0.40)

# Accepted equilibrium connectance range
C_eq_range = (0.03, 0.35)

Z = 10.0

tmax = 10_000
tmin = 1_000

extinction_threshold = 1e-12
Re_ode_tolerance = 1e-6

checkpoint_every_accepts = 1
checkpoint_every_attempts = 500

# Initial-richness proposal windows for each intended
# realised-richness bin.
S_proposal_override = Dict(
    "S6_10"  => (10, 35),
    "S11_15" => (14, 50),
    "S16_20" => (20, 70),
    "S21_25" => (28, 90),
    "S26_30" => (40, 110),
    "S31_35" => (52, 130),
    "S36_40" => (65, 150)
)

out_dir = "data/model_fws/classic_response_h2_sim_7bins_70tasks"
mkpath(out_dir)

out_file = joinpath(
    out_dir,
    "model_food_webs_task_$(task_tag).jld2"
)

checkpoint_file = joinpath(
    out_dir,
    "model_food_webs_task_$(task_tag)_checkpoint.jld2"
)

# ============================================================
# Helper functions
# ============================================================

function add_fail!(fail_counts, key)
    fail_counts[key] = get(fail_counts, key, 0) + 1
    return nothing
end

function empty_bin_counts(S_bins)
    return Dict(b.label => 0 for b in S_bins)
end

function allocate_task_replicates(
    S_bins,
    n_rep,
    task_id,
    n_tasks
)
    return Dict(
        b.label => collect(task_id:n_tasks:n_rep)
        for b in S_bins
    )
end

function task_target_counts(task_replicates)
    return Dict(
        label => length(reps)
        for (label, reps) in task_replicates
    )
end

function match_S_bin(
    S_eq,
    S_bins,
    counts,
    task_targets
)
    for b in S_bins
        in_bin = b.min <= S_eq <= b.max
        bin_not_full = counts[b.label] < task_targets[b.label]

        if in_bin && bin_not_full
            return b
        end
    end

    return nothing
end

function proposal_window(
    S_bin,
    S_range,
    overrides
)
    if haskey(overrides, S_bin.label)
        S_low, S_high = overrides[S_bin.label]
    else
        S_low = floor(Int, S_bin.target / 0.95)
        S_high = ceil(Int, S_bin.target / 0.40)
    end

    S_low = max(first(S_range), S_low)
    S_high = min(last(S_range), S_high)

    return S_low, S_high
end

function propose_initial_web(
    S_bins,
    counts,
    task_targets,
    S_range,
    C_range,
    overrides
)
    open_bins = [
        b for b in S_bins
        if counts[b.label] < task_targets[b.label]
    ]

    isempty(open_bins) && return nothing

    missing = [
        task_targets[b.label] - counts[b.label]
        for b in open_bins
    ]

    proposed_bin = sample(
        open_bins,
        Weights(missing)
    )

    S_low, S_high = proposal_window(
        proposed_bin,
        S_range,
        overrides
    )

    S_low > S_high && return nothing

    S = rand(S_low:S_high)

    C_low = max(
        C_range[1],
        (S - 1) / S^2
    )

    C_low >= C_range[2] && return nothing

    C = rand(Uniform(C_low, C_range[2]))

    return S, C, proposed_bin.label
end

function classify_error(err)
    msg = lowercase(sprint(showerror, err))

    occursin("could not generate adequate network", msg) &&
        return "niche_generation_failed"

    occursin("connectance", msg) &&
        return "connectance"

    (
        occursin("nan", msg) ||
        occursin("inf", msg)
    ) && return "nan_or_inf"

    occursin("dt", msg) &&
        return "solver_dt_failure"

    return "other"
end

function print_progress(
    task_id,
    attempt,
    attempt_start,
    counts,
    targets,
    fail_counts
)
    println(
        "Task $task_id | attempt $attempt | ",
        "this run $(attempt - attempt_start)"
    )

    println(
        "Accepted: ",
        sum(values(counts)),
        " / ",
        sum(values(targets))
    )

    println("Counts: ", counts)
    println("Failures: ", fail_counts)

    flush(stdout)
end

function save_checkpoint(
    file;
    task_id,
    n_tasks,
    attempt,
    valid_food_webs,
    S_bins,
    task_replicates,
    task_targets,
    target_counts,
    fail_counts,
    n_rep_per_Sbin,
    S_range,
    S_proposal_override,
    hill_exponent,
    i0_intra
)
    jldsave(
        file;
        task_id,
        n_tasks,
        attempt,
        valid_food_webs,
        S_bins,
        task_replicates,
        task_targets,
        target_counts,
        fail_counts,
        n_rep_per_Sbin,
        S_range,
        S_proposal_override,
        hill_exponent,
        i0_intra
    )

    println("Checkpoint saved at attempt $attempt")
    flush(stdout)
end

function save_final(
    file;
    task_id,
    n_tasks,
    attempt,
    valid_food_webs,
    S_bins,
    task_replicates,
    task_targets,
    target_counts,
    fail_counts,
    n_rep_per_Sbin,
    S_range,
    S_proposal_override,
    C_eq_range,
    Re_ode_tolerance,
    hill_exponent,
    i0_intra
)
    jldsave(
        file;
        task_id,
        n_tasks,
        attempt,
        valid_food_webs,
        S_bins,
        task_replicates,
        task_targets,
        target_counts,
        fail_counts,
        n_rep_per_Sbin,
        S_range,
        S_proposal_override,
        C_eq_range,
        Re_ode_tolerance,
        hill_exponent,
        i0_intra
    )
end

# ============================================================
# Task targets
# ============================================================

task_replicates = allocate_task_replicates(
    S_bins,
    n_rep_per_Sbin,
    task_id,
    n_tasks
)

task_targets = task_target_counts(task_replicates)

n_valid_target_task = sum(values(task_targets))

println("Starting task $task_id / $n_tasks")
println("Task replicates: ", task_replicates)
println("Task target: $n_valid_target_task webs")
println("Proposal windows: ", S_proposal_override)
println("Hill exponent: ", hill_exponent)
println("Intraspecific interference: ", i0_intra)

flush(stdout)

# ============================================================
# Resume or start
# ============================================================

if isfile(checkpoint_file)
    checkpoint = load(checkpoint_file)

    checkpoint["task_id"] == task_id ||
        error("Checkpoint belongs to another task")

    checkpoint["n_tasks"] == n_tasks ||
        error("Checkpoint used a different n_tasks")

    checkpoint["n_rep_per_Sbin"] == n_rep_per_Sbin ||
        error("Checkpoint used a different n_rep_per_Sbin")

    checkpoint_labels = Set(
        b.label for b in checkpoint["S_bins"]
    )

    current_labels = Set(
        b.label for b in S_bins
    )

    checkpoint_labels == current_labels ||
        error("Checkpoint used a different set of S bins")

    checkpoint["S_range"] == S_range ||
        error("Checkpoint used a different S_range")

    checkpoint["S_proposal_override"] == S_proposal_override ||
        error("Checkpoint used different proposal windows")

    checkpoint["hill_exponent"] == hill_exponent ||
        error("Checkpoint used a different hill_exponent")

    checkpoint["i0_intra"] == i0_intra ||
        error("Checkpoint used a different i0_intra")

    valid_food_webs = checkpoint["valid_food_webs"]
    target_counts = checkpoint["target_counts"]
    fail_counts = checkpoint["fail_counts"]
    attempt = checkpoint["attempt"]

    println(
        "Resuming from attempt $attempt with ",
        "$(sum(values(target_counts))) accepted webs"
    )
else
    valid_food_webs = Any[]
    target_counts = empty_bin_counts(S_bins)
    fail_counts = Dict{String,Int}()
    attempt = 0

    println("Starting without checkpoint")
end

flush(stdout)

attempt_start = attempt
accepted_at_start = sum(values(target_counts))

# ============================================================
# Main simulation
# ============================================================

while (
    sum(values(target_counts)) < n_valid_target_task &&
    attempt - attempt_start < max_attempts_this_run
)
    global attempt += 1

    Random.seed!(
        base_seed +
        task_id * 10_000_000 +
        attempt
    )

    try
        proposed = propose_initial_web(
            S_bins,
            target_counts,
            task_targets,
            S_range,
            C_range,
            S_proposal_override
        )

        if isnothing(proposed)
            add_fail!(fail_counts, "no_proposal")
            continue
        end

        S, C, proposed_S_bin = proposed

        # ----------------------------------------------------
        # Initial niche-model food web
        # ----------------------------------------------------

        fw = Foodweb(
            :niche;
            S=S,
            C=C,
            reject_cycles=true,
            tol_C=0.1
        )

        A = get_A(fw)

        if sum(A .> 0) == 0
            add_fail!(
                fail_counts,
                "empty_initial_web"
            )

            continue
        end

        # ----------------------------------------------------
        # Burn-in simulation
        # ----------------------------------------------------

        model_burn = default_model(
            fw,
            ClassicResponse(
                ;
                h=hill_exponent,
                i0=i0_intra
            ),
            BodyMass(; Z=Z),
            Efficiency(:Miele2019)
        )

        B0 = rand(S) .+ 0.1

        sol_burn = simulate(
            model_burn,
            B0,
            tmax;
            callback=CallbackSet(
                extinction_callback(
                    model_burn,
                    extinction_threshold;
                    verbose=false
                ),
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

        if !all(isfinite, B_burn)
            add_fail!(
                fail_counts,
                "burn_nonfinite_B"
            )

            continue
        end

        A_burn0 = get_A(model_burn)

        is_producer_burn0 = vec(
            sum(
                A_burn0 .> 0;
                dims=2
            ) .== 0
        )

        alive_burn = findall(
            B_burn .> extinction_threshold
        )

        if isempty(alive_burn)
            add_fail!(
                fail_counts,
                "no_alive_after_burn"
            )

            continue
        end

        keep_burn = prune_realised_species(
            A_burn0,
            alive_burn,
            is_producer_burn0
        )

        if isempty(keep_burn)
            add_fail!(
                fail_counts,
                "no_realised_species_after_burn"
            )

            continue
        end

        A_burn = A_burn0[
            keep_burn,
            keep_burn
        ]

        if sum(A_burn .> 0) == 0
            add_fail!(
                fail_counts,
                "empty_realised_web_after_burn"
            )

            continue
        end

        B0_finalrun = B_burn[keep_burn]
        fw_burn = Foodweb(A_burn)

        attack_rate_burn =
            Matrix(model_burn.attack_rate)[
                keep_burn,
                keep_burn
            ]

        handling_time_burn =
            Matrix(model_burn.handling_time)[
                keep_burn,
                keep_burn
            ]

        preference_burn =
            Matrix(model_burn.consumers_preferences)[
                keep_burn,
                keep_burn
            ]

        efficiency_burn =
            Matrix(model_burn.efficiency)[
                keep_burn,
                keep_burn
            ]

        body_mass_burn =
            collect(model_burn.body_mass)[
                keep_burn
            ]

        # ----------------------------------------------------
        # Final equilibrium simulation
        # ----------------------------------------------------

        model = default_model(
            fw_burn,
            ClassicResponse(
                ;
                h=hill_exponent,
                a=attack_rate_burn,
                h_t=handling_time_burn,
                w=preference_burn,
                i0=i0_intra
            ),
            BodyMass(body_mass_burn),
            Efficiency(efficiency_burn)
        )

        sol = simulate(
            model,
            B0_finalrun,
            tmax;
            callback=CallbackSet(
                extinction_callback(
                    model,
                    extinction_threshold;
                    verbose=false
                ),
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

        if !all(isfinite, B_final)
            add_fail!(
                fail_counts,
                "final_nonfinite_B"
            )

            continue
        end

        if !all(B_final .> extinction_threshold)
            add_fail!(
                fail_counts,
                "extinction_after_final_run"
            )

            continue
        end

        # ----------------------------------------------------
        # Extract final-model parameters
        # ----------------------------------------------------

        M = collect(model.body_mass)

        a = Matrix(model.attack_rate)
        ht = Matrix(model.handling_time)
        w = Matrix(model.consumers_preferences)
        e = Matrix(model.efficiency)

        r = collect(model.growth_rate)
        K = collect(model.carrying_capacity)
        x = collect(model.metabolism)

        # ----------------------------------------------------
        # Equilibrium summary and Jacobians
        # ----------------------------------------------------

        eq = equilibrium_foodweb_summary_classic(
            fw_burn,
            B_final,
            M,
            a,
            ht,
            w,
            e,
            hill_exponent,
            i0_intra,
            r,
            K,
            x;
            extinction_threshold=extinction_threshold
        )

        if isnothing(eq)
            add_fail!(
                fail_counts,
                "no_equilibrium_summary"
            )

            continue
        end

        if eq.Re_ode_eq > Re_ode_tolerance
            add_fail!(
                fail_counts,
                "unstable_ode"
            )

            continue
        end

        S_bin = match_S_bin(
            eq.S_eq,
            S_bins,
            target_counts,
            task_targets
        )

        if isnothing(S_bin)
            add_fail!(
                fail_counts,
                "wrong_S_eq_bin_or_bin_full"
            )

            continue
        end

        if !(
            C_eq_range[1] <=
            eq.C_eq <=
            C_eq_range[2]
        )
            add_fail!(
                fail_counts,
                "wrong_C_eq"
            )

            continue
        end

        # ----------------------------------------------------
        # Accept food web
        # ----------------------------------------------------

        next_position =
            target_counts[S_bin.label] + 1

        rep = task_replicates[
            S_bin.label
        ][next_position]

        bin_index = findfirst(
            b -> b.label == S_bin.label,
            S_bins
        )

        valid_id =
            (bin_index - 1) *
            n_rep_per_Sbin +
            rep

        target_counts[S_bin.label] += 1

        accepted_total =
            sum(values(target_counts))

        push!(
            valid_food_webs,
            (
                valid_id=valid_id,
                task_id=task_id,

                S_label=S_bin.label,
                S_min=S_bin.min,
                S_max=S_bin.max,
                S_target=S_bin.target,
                rep=rep,

                initial_S=S,
                initial_C=C,
                proposed_S_bin=proposed_S_bin,

                initial_A=A,
                initial_L=sum(A .> 0),

                functional_response="ClassicResponse",
                hill_exponent=hill_exponent,
                i0_intra=i0_intra,

                body_mass=M,
                attack_rate=a,
                handling_time=ht,
                preference=w,
                efficiency=e,

                r=r,
                K=K,
                x=x,

                B_final=B_final,

                alive_persistence=persistence(
                    B_final;
                    threshold=extinction_threshold
                ),

                final_persistence=eq.S_eq / S,

                keep_species=eq.keep_species,

                S_eq=eq.S_eq,
                L_eq=eq.L_eq,
                C_eq=eq.C_eq,
                maxTL_eq=eq.maxTL_eq,

                A_eq=eq.A_eq,
                B_eq=eq.B_eq,
                TL_eq=eq.TL_eq,
                is_producer_eq=eq.is_producer_eq,

                M_eq=eq.M_eq,

                a_eq=eq.a_eq,
                ht_eq=eq.ht_eq,
                w_eq=eq.w_eq,
                e_eq=eq.e_eq,

                hill_exponent_eq=eq.h_eq,
                i0_eq=eq.i0_eq,

                r_eq=eq.r_eq,
                K_eq=eq.K_eq,
                x_eq=eq.x_eq,

                denominator_eq=eq.denominator_eq,
                saturation_load_eq=eq.saturation_load_eq,
                saturation_index_eq=eq.saturation_index_eq,

                feeding_rate_eq=eq.feeding_rate_eq,
                consumption_flux_eq=eq.consumption_flux_eq,
                assimilated_flux_eq=eq.assimilated_flux_eq,

                J_emp=eq.J_emp_eq,
                J_ode=eq.J_ode_eq,
                J_ode_diag=eq.J_ode_diag,

                Re_emp=eq.Re_emp_eq,
                Re_ode=eq.Re_ode_eq,

                positive_diag_idx=eq.positive_diag_idx,
                positive_diag_values=eq.positive_diag_values,

                totalB_final=sum(B_final),
                totalB_eq=sum(eq.B_eq)
            )
        )

        println(
            "Accepted $accepted_total / ",
            "$n_valid_target_task | ",
            "bin=$(S_bin.label) | ",
            "rep=$rep | ",
            "valid_id=$valid_id | ",
            "proposed_bin=$proposed_S_bin | ",
            "initial_S=$S | ",
            "S_eq=$(eq.S_eq) | ",
            "persistence=$(round(eq.S_eq / S, digits=3)) | ",
            "C_eq=$(round(eq.C_eq, digits=3)) | ",
            "Re_ode=$(round(eq.Re_ode_eq, digits=5)) | ",
            "Re_emp=$(round(eq.Re_emp_eq, digits=5))"
        )

        flush(stdout)

        if accepted_total %
           checkpoint_every_accepts == 0

            save_checkpoint(
                checkpoint_file;
                task_id=task_id,
                n_tasks=n_tasks,
                attempt=attempt,
                valid_food_webs=valid_food_webs,
                S_bins=S_bins,
                task_replicates=task_replicates,
                task_targets=task_targets,
                target_counts=target_counts,
                fail_counts=fail_counts,
                n_rep_per_Sbin=n_rep_per_Sbin,
                S_range=S_range,
                S_proposal_override=S_proposal_override,
                hill_exponent=hill_exponent,
                i0_intra=i0_intra
            )
        end

    catch err
        add_fail!(
            fail_counts,
            classify_error(err)
        )
    end

    if attempt %
       checkpoint_every_attempts == 0

        save_checkpoint(
            checkpoint_file;
            task_id=task_id,
            n_tasks=n_tasks,
            attempt=attempt,
            valid_food_webs=valid_food_webs,
            S_bins=S_bins,
            task_replicates=task_replicates,
            task_targets=task_targets,
            target_counts=target_counts,
            fail_counts=fail_counts,
            n_rep_per_Sbin=n_rep_per_Sbin,
            S_range=S_range,
            S_proposal_override=S_proposal_override,
            hill_exponent=hill_exponent,
            i0_intra=i0_intra
        )

        print_progress(
            task_id,
            attempt,
            attempt_start,
            target_counts,
            task_targets,
            fail_counts
        )
    end
end

# ============================================================
# Final save
# ============================================================

save_checkpoint(
    checkpoint_file;
    task_id=task_id,
    n_tasks=n_tasks,
    attempt=attempt,
    valid_food_webs=valid_food_webs,
    S_bins=S_bins,
    task_replicates=task_replicates,
    task_targets=task_targets,
    target_counts=target_counts,
    fail_counts=fail_counts,
    n_rep_per_Sbin=n_rep_per_Sbin,
    S_range=S_range,
    S_proposal_override=S_proposal_override,
    hill_exponent=hill_exponent,
    i0_intra=i0_intra
)

save_final(
    out_file;
    task_id=task_id,
    n_tasks=n_tasks,
    attempt=attempt,
    valid_food_webs=valid_food_webs,
    S_bins=S_bins,
    task_replicates=task_replicates,
    task_targets=task_targets,
    target_counts=target_counts,
    fail_counts=fail_counts,
    n_rep_per_Sbin=n_rep_per_Sbin,
    S_range=S_range,
    S_proposal_override=S_proposal_override,
    C_eq_range=C_eq_range,
    Re_ode_tolerance=Re_ode_tolerance,
    hill_exponent=hill_exponent,
    i0_intra=i0_intra
)

n_accepted = sum(values(target_counts))

println("================================================")
println("Task $task_id finished")
println("Complete: ", n_accepted == n_valid_target_task)

println(
    "Accepted this run: ",
    n_accepted - accepted_at_start
)

println(
    "Total accepted: ",
    "$n_accepted / $n_valid_target_task"
)

println("Counts: ", target_counts)
println("Failures: ", fail_counts)
println("Cumulative attempts: ", attempt)
println("Saved: ", out_file)
println("Checkpoint: ", checkpoint_file)
println("================================================")

flush(stdout)