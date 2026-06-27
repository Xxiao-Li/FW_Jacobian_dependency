const ENDI = EcologicalNetworksDynamics.Internals

get_A(fw) = Matrix(fw.A)

function random_linear_response_params(fw; alpha_range = (1.0, 10.0), α_range = nothing)
    if !isnothing(α_range)
        alpha_range = α_range
    end

    A = get_A(fw)
    S = size(A, 1)

    alpha = zeros(Float64, S)
    omega = zeros(Float64, S, S)

    alpha_dist = Uniform(alpha_range[1], alpha_range[2])

    for consumer in 1:S
        prey_list = findall(A[consumer, :] .> 0)

        if !isempty(prey_list)
            alpha[consumer] = rand(alpha_dist)

            raw_omega = rand(length(prey_list))
            omega[consumer, prey_list] .= raw_omega ./ sum(raw_omega)
        end
    end

    return alpha, sparse(omega)
end

function persistence(B; threshold = 1e-8)
    mean(B .> threshold)
end

function re_lmax(J)
    maximum(real.(eigvals(J)))
end

function prune_realised_species(A, alive, is_producer_original)
    A_alive = A[alive, alive]
    is_producer_alive = is_producer_original[alive]

    keep_local = trues(length(alive))

    while true
        local_idx = findall(keep_local)
        isempty(local_idx) && return Int[]

        A_sub = A_alive[local_idx, local_idx]
        is_producer_sub = is_producer_alive[local_idx]

        has_prey = vec(sum(A_sub .> 0, dims = 2)) .> 0
        has_predator = vec(sum(A_sub .> 0, dims = 1)) .> 0

        # Basal producers must support at least one consumer.
        # Consumers must retain at least one prey item.
        retain_sub =
            (is_producer_sub .& has_predator) .|
            (.!is_producer_sub .& has_prey)

        new_keep_local = falses(length(alive))
        new_keep_local[local_idx[retain_sub]] .= true

        new_keep_local == keep_local && break
        keep_local = new_keep_local
    end

    return alive[findall(keep_local)]
end

function jacobian_empirical(fw, B, alpha, omega, e)
    A = get_A(fw)
    S = size(A, 1)

    J = zeros(Float64, S, S)

    for consumer in 1:S
        prey_list = findall(A[consumer, :] .> 0)

        for resource in prey_list
            c = alpha[consumer] * omega[consumer, resource]

            J[resource, consumer] += -c * B[resource]
            J[consumer, resource] += e[consumer, resource] * c * B[consumer]
        end
    end

    J[diagind(J)] .= 0.0

    return J
end

function jacobian_ode(fw, B, alpha, omega, e, r, K, x, is_producer)
    A = get_A(fw)
    S = size(A, 1)

    J = zeros(Float64, S, S)

    for consumer in 1:S
        prey_list = findall(A[consumer, :] .> 0)

        for resource in prey_list
            c = alpha[consumer] * omega[consumer, resource]

            # Resource equation:
            # dB_resource/dt includes -B_consumer * c * B_resource.
            J[resource, consumer] += -c * B[resource]
            J[resource, resource] += -c * B[consumer]

            # Consumer equation:
            # dB_consumer/dt includes e * B_consumer * c * B_resource.
            J[consumer, resource] += e[consumer, resource] * c * B[consumer]
            J[consumer, consumer] += e[consumer, resource] * c * B[resource]
        end
    end

    for i in 1:S
        if is_producer[i]
            J[i, i] += r[i] * (1 - 2 * B[i] / K[i])
        end

        J[i, i] += -x[i]
    end

    return J
end

function equilibrium_foodweb_summary(fw, B_final, alpha, omega, e, r, K, x;
    extinction_threshold = 1e-8
)
    A = get_A(fw)

    # Producer identity must come from the original ODE model, not from the
    # pruned realised web.
    is_producer_original = vec(sum(A .> 0, dims = 2) .== 0)

    alive = findall(B_final .> extinction_threshold)
    isempty(alive) && return nothing

    keep = prune_realised_species(A, alive, is_producer_original)
    isempty(keep) && return nothing

    A_eq = A[keep, keep]
    L_eq = sum(A_eq .> 0)
    L_eq == 0 && return nothing

    B_eq = B_final[keep]

    alpha_eq = alpha[keep]
    omega_eq = Matrix(omega)[keep, keep]

    e_eq = e[keep, keep]
    r_eq = r[keep]
    K_eq = K[keep]
    x_eq = x[keep]
    is_producer_eq = is_producer_original[keep]

    S_eq = length(keep)
    C_eq = L_eq / S_eq^2

    TL_eq = ENDI.trophic_levels(A_eq)
    maxTL_eq = maximum(TL_eq)

    fw_eq = Foodweb(A_eq)

    J_emp_eq = jacobian_empirical(
        fw_eq,
        B_eq,
        alpha_eq,
        omega_eq,
        e_eq
    )

    J_ode_eq = jacobian_ode(
        fw_eq,
        B_eq,
        alpha_eq,
        omega_eq,
        e_eq,
        r_eq,
        K_eq,
        x_eq,
        is_producer_eq
    )

    Re_emp_eq = re_lmax(J_emp_eq)
    Re_ode_eq = re_lmax(J_ode_eq)

    return (
        keep_species = keep,
        S_eq = S_eq,
        L_eq = L_eq,
        C_eq = C_eq,
        maxTL_eq = maxTL_eq,
        A_eq = A_eq,
        B_eq = B_eq,
        TL_eq = TL_eq,
        is_producer_eq = is_producer_eq,
        alpha_eq = alpha_eq,
        omega_eq = sparse(omega_eq),
        α_eq = alpha_eq,
        ω_eq = sparse(omega_eq),
        e_eq = e_eq,
        r_eq = r_eq,
        K_eq = K_eq,
        x_eq = x_eq,
        J_emp_eq = J_emp_eq,
        J_ode_eq = J_ode_eq,
        Re_emp_eq = Re_emp_eq,
        Re_ode_eq = Re_ode_eq,
        positive_diag_idx = findall(diag(J_ode_eq) .> 0),
        positive_diag_values = diag(J_ode_eq)[findall(diag(J_ode_eq) .> 0)]
    )
end
