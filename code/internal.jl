get_A(fw) = Matrix(fw.A)

function random_linear_response_params(fw; α_range = (1.0, 10.0))
    A = get_A(fw)
    S = size(A, 1)

    α = zeros(Float64, S)
    ω = zeros(Float64, S, S)

    α_dist = Uniform(α_range[1], α_range[2])

    for consumer in 1:S
        prey_list = findall(A[consumer, :] .> 0)

        # Only consumers receive feeding parameters
        if !isempty(prey_list)
            # αᵢ ranges from 1 to 10
            α[consumer] = rand(α_dist)

            # 0 ≤ ωᵢⱼ ≤ 1 and sum_j ωᵢⱼ = 1 for the same consumer
            raw_ω = rand(length(prey_list))
            ω[consumer, prey_list] .= raw_ω ./ sum(raw_ω)
        end
    end

    return α, sparse(ω)
end

function persistence(B; threshold = 1e-8)
    mean(B .> threshold)
end

function re_lmax(J)
    maximum(real.(eigvals(J)))
end


function jacobian_empirical(fw, B, α, ω, e)
    A = get_A(fw)
    S = size(A, 1)

    J = zeros(Float64, S, S)

    for consumer in 1:S
        prey_list = findall(A[consumer, :] .> 0)

        for resource in prey_list
            c = α[consumer] * ω[consumer, resource]

            # Effect of consumer on resource
            # J[resource, consumer] = -cᵢⱼ B_resource
            J[resource, consumer] += -c * B[resource]

            # Effect of resource on consumer
            # J[consumer, resource] = eᵢⱼ cᵢⱼ B_consumer
            J[consumer, resource] += e[consumer, resource] * c * B[consumer]
        end
    end

    # Empirical Jacobian assumes no diagonal self-effects
    for i in 1:S
        J[i, i] = 0.0
    end

    return J
end

function jacobian_ode(fw, B, α, ω, e, r, K, x)
    A = get_A(fw)
    S = size(A, 1)

    J = zeros(Float64, S, S)

    for consumer in 1:S
        prey_list = findall(A[consumer, :] .> 0)

        for resource in prey_list
            c = α[consumer] * ω[consumer, resource]

            # ------------------------------------------------
            # Resource equation:
            # dB_resource/dt includes -B_consumer * c * B_resource
            # ------------------------------------------------
            J[resource, consumer] += -c * B[resource]
            J[resource, resource] += -c * B[consumer]

            # ------------------------------------------------
            # Consumer equation:
            # dB_consumer/dt includes +e * B_consumer * c * B_resource
            # ------------------------------------------------
            J[consumer, resource] += e[consumer, resource] * c * B[consumer]
            J[consumer, consumer] += e[consumer, resource] * c * B[resource]
        end
    end

    # Producer logistic growth:
    # rᵢBᵢ(1 - Bᵢ/Kᵢ)
    for i in 1:S
        is_consumer = any(A[i, :] .> 0)

        if !is_consumer
            J[i, i] += r[i] * (1 - 2 * B[i] / K[i])
        end
    end

    # Metabolism and mortality:
    # -xᵢBᵢ - dᵢBᵢ
    for i in 1:S
        J[i, i] += -x[i]
    end

    return J
end


function equilibrium_foodweb_summary(fw, B_final, α, ω, e, r, K, x;
    extinction_threshold = 1e-8
)
    A = get_A(fw)
    S = size(A, 1)

    # --------------------------------------------------------
    # 1. Keep alive species first
    # --------------------------------------------------------

    alive = findall(B_final .> extinction_threshold)

    if length(alive) == 0
        return nothing
    end

    A_alive = A[alive, alive]

    # --------------------------------------------------------
    # 2. Remove disconnected species within the alive subweb
    # A species is connected if it has at least one incoming or outgoing link
    # --------------------------------------------------------

    row_degree = vec(sum(A_alive .> 0, dims = 2))   # consumer links
    col_degree = vec(sum(A_alive .> 0, dims = 1))   # prey links

    connected_local = findall((row_degree .+ col_degree) .> 0)

    if length(connected_local) == 0
        return nothing
    end

    keep = alive[connected_local]

    # --------------------------------------------------------
    # 3. Subset equilibrium food web and parameters
    # --------------------------------------------------------

    A_eq = A[keep, keep]

    B_eq = B_final[keep]

    α_eq = α[keep]
    ω_eq = Matrix(ω)[keep, keep]

    e_eq = e[keep, keep]
    r_eq = r[keep]
    K_eq = K[keep]
    x_eq = x[keep]

    S_eq = length(keep)
    L_eq = sum(A_eq .> 0)
    C_eq = L_eq / S_eq^2

    # --------------------------------------------------------
    # 4. Build a new food web object for the realised subweb
    # --------------------------------------------------------

    fw_eq = Foodweb(A_eq)

    # --------------------------------------------------------
    # 5. Calculate Jacobian matrices for the realised subweb
    # --------------------------------------------------------

    J_emp_eq = jacobian_empirical(
        fw_eq,
        B_eq,
        α_eq,
        ω_eq,
        e_eq
    )

    J_ode_eq = jacobian_ode(
        fw_eq,
        B_eq,
        α_eq,
        ω_eq,
        e_eq,
        r_eq,
        K_eq,
        x_eq
    )

    Re_emp_eq = re_lmax(J_emp_eq)
    Re_ode_eq = re_lmax(J_ode_eq)

    return (
        keep_species = keep,
        S_eq = S_eq,
        L_eq = L_eq,
        C_eq = C_eq,
        A_eq = A_eq,
        B_eq = B_eq,
        α_eq = α_eq,
        ω_eq = sparse(ω_eq),
        e_eq = e_eq,
        r_eq = r_eq,
        K_eq = K_eq,
        x_eq = x_eq,
        J_emp_eq = J_emp_eq,
        J_ode_eq = J_ode_eq,
        Re_emp_eq = Re_emp_eq,
        Re_ode_eq = Re_ode_eq
    )
end