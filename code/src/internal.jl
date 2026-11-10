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

    # Producer identity must come from the original ODE model.
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




function classic_equilibrium_quantities(
    A, B, M, a, ht, w, e, h, i0
)
    S = length(B)
    i0 = i0 isa Number ? fill(i0, S) : collect(i0)

    denominator = ones(S)
    feeding_rate = zeros(S, S)
    consumption_flux = zeros(S, S)
    assimilated_flux = zeros(S, S)

    for consumer in 1:S
        prey = findall(A[consumer, :] .> 0)
        isempty(prey) && continue

        denominator[consumer] += i0[consumer] * B[consumer]

        for resource in prey
            denominator[consumer] +=
                w[consumer, resource] *
                a[consumer, resource] *
                ht[consumer, resource] *
                B[resource]^h
        end

        for resource in prey
            feeding_rate[consumer, resource] =
                w[consumer, resource] *
                a[consumer, resource] *
                B[resource]^h /
                (
                    M[consumer] *
                    denominator[consumer]
                )

            consumption_flux[consumer, resource] =
                B[consumer] *
                feeding_rate[consumer, resource]

            assimilated_flux[consumer, resource] =
                e[consumer, resource] *
                consumption_flux[consumer, resource]
        end
    end

    saturation_load = denominator .- 1
    saturation_index = saturation_load ./ denominator

    return (
        denominator=denominator,
        saturation_load=saturation_load,
        saturation_index=saturation_index,
        feeding_rate=feeding_rate,
        consumption_flux=consumption_flux,
        assimilated_flux=assimilated_flux
    )
end


function classic_response_derivatives(
    A, B, M, a, ht, w, h, i0, denominator
)
    S = length(B)
    i0 = i0 isa Number ? fill(i0, S) : collect(i0)

    dF = zeros(S, S, S)

    for consumer in 1:S
        prey = findall(A[consumer, :] .> 0)
        isempty(prey) && continue

        D = denominator[consumer]

        for resource in prey
            q =
                w[consumer, resource] *
                a[consumer, resource]

            numerator =
                q *
                B[resource]^h /
                M[consumer]

            for species in 1:S
                dnumerator =
                    species == resource ?
                    q * h * B[resource]^(h - 1) /
                    M[consumer] :
                    0.0

                ddenominator =
                    species == consumer ?
                    i0[consumer] :
                    0.0

                if A[consumer, species] > 0
                    ddenominator +=
                        w[consumer, species] *
                        a[consumer, species] *
                        ht[consumer, species] *
                        h *
                        B[species]^(h - 1)
                end

                dF[consumer, resource, species] =
                    dnumerator / D -
                    numerator * ddenominator / D^2
            end
        end
    end

    return dF
end


function jacobian_ode_classic(
    A, B, M, a, ht, w, e, h, i0,
    r, K, x, is_producer
)
    S = length(B)

    q = classic_equilibrium_quantities(
        A, B, M, a, ht, w, e, h, i0
    )

    dF = classic_response_derivatives(
        A, B, M, a, ht, w, h, i0,
        q.denominator
    )

    J = zeros(S, S)

    for consumer in 1:S
        prey = findall(A[consumer, :] .> 0)

        for resource in prey
            F = q.feeding_rate[consumer, resource]

            for species in 1:S
                dflux =
                    B[consumer] *
                    dF[consumer, resource, species]

                if species == consumer
                    dflux += F
                end

                J[resource, species] -= dflux

                J[consumer, species] +=
                    e[consumer, resource] *
                    dflux
            end
        end
    end

    for i in 1:S
        if is_producer[i]
            J[i, i] +=
                r[i] *
                (1 - 2 * B[i] / K[i])
        end

        J[i, i] -= x[i]
    end

    return J, q
end


function equilibrium_foodweb_summary_classic(
    fw,
    B_final,
    M,
    a,
    ht,
    w,
    e,
    h,
    i0,
    r,
    K,
    x;
    extinction_threshold=1e-8
)
    A = get_A(fw)

    is_producer = vec(
        sum(A .> 0, dims=2) .== 0
    )

    alive = findall(
        B_final .> extinction_threshold
    )

    isempty(alive) && return nothing

    keep = prune_realised_species(
        A,
        alive,
        is_producer
    )

    isempty(keep) && return nothing

    A_eq = A[keep, keep]
    B_eq = B_final[keep]

    L_eq = sum(A_eq .> 0)
    L_eq == 0 && return nothing

    M_eq = M[keep]
    a_eq = Matrix(a)[keep, keep]
    ht_eq = Matrix(ht)[keep, keep]
    w_eq = Matrix(w)[keep, keep]
    e_eq = Matrix(e)[keep, keep]

    r_eq = r[keep]
    K_eq = K[keep]
    x_eq = x[keep]

    i0_full =
        i0 isa Number ?
        fill(i0, length(B_final)) :
        collect(i0)

    i0_eq = i0_full[keep]
    is_producer_eq = is_producer[keep]

    J_ode_eq, q = jacobian_ode_classic(
        A_eq,
        B_eq,
        M_eq,
        a_eq,
        ht_eq,
        w_eq,
        e_eq,
        h,
        i0_eq,
        r_eq,
        K_eq,
        x_eq,
        is_producer_eq
    )

    J_emp_eq = copy(J_ode_eq)
    J_emp_eq[diagind(J_emp_eq)] .= 0

    S_eq = length(keep)
    C_eq = L_eq / S_eq^2

    TL_eq = ENDI.trophic_levels(A_eq)
    J_diag = diag(J_ode_eq)
    positive_diag_idx = findall(J_diag .> 0)

    return (
        keep_species=keep,

        S_eq=S_eq,
        L_eq=L_eq,
        C_eq=C_eq,
        maxTL_eq=maximum(TL_eq),

        A_eq=A_eq,
        B_eq=B_eq,
        TL_eq=TL_eq,
        is_producer_eq=is_producer_eq,

        M_eq=M_eq,
        a_eq=sparse(a_eq),
        ht_eq=sparse(ht_eq),
        w_eq=sparse(w_eq),
        e_eq=sparse(e_eq),

        h_eq=h,
        i0_eq=i0_eq,

        r_eq=r_eq,
        K_eq=K_eq,
        x_eq=x_eq,

        denominator_eq=q.denominator,
        saturation_load_eq=q.saturation_load,
        saturation_index_eq=q.saturation_index,

        feeding_rate_eq=sparse(q.feeding_rate),
        consumption_flux_eq=sparse(q.consumption_flux),
        assimilated_flux_eq=sparse(q.assimilated_flux),

        J_emp_eq=J_emp_eq,
        J_ode_eq=J_ode_eq,
        J_ode_diag=J_diag,

        Re_emp_eq=re_lmax(J_emp_eq),
        Re_ode_eq=re_lmax(J_ode_eq),

        positive_diag_idx=positive_diag_idx,
        positive_diag_values=
            J_diag[positive_diag_idx]
    )
end