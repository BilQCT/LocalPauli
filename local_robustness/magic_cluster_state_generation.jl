# Historical E(G) T^U |+> resource convention; NOT the manuscript robustness target.
# Require an explicit opt-in to prevent accidentally generating publication data.
if abspath(PROGRAM_FILE) == @__FILE__ && !("--legacy" in ARGS)
    error("Historical state convention. Use the bicolored manuscript pipeline, or pass --legacy explicitly.")
end
using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using JSON3, HDF5, LinearAlgebra

data_folder = joinpath(base_dir, "data")
libs_folder = joinpath(base_dir, "libs")

push!(LOAD_PATH, libs_folder)

function all_dit_strings(base::Int, n::Int)
    return [digits(i, base=base, pad=n) for i in 0:(base^n - 1)]
end

include(joinpath(libs_folder, "pauli.jl"))

# ------------------------------------------------------------------------------
# Pauli Basis Actions (from your existing code)
# ------------------------------------------------------------------------------

function hadamard_action_on_bsf(q, a)
    n = div(length(a), 2)
    a_prime = copy(a)
    a_prime[q] = a[n+q]
    a_prime[n+q] = a[q]
    r = (a[q] * a[n+q]) % 2
    return a_prime, r
end

function phase_action_on_bsf(q, a)
    n = div(length(a), 2)
    a_prime = copy(a)
    a_prime[n+q] = (a[q] + a[n+q]) % 2
    r = (a[q] * a[n+q]) % 2
    return a_prime, r
end

function cx_action_on_bsf(control, target, a)
    n = div(length(a), 2)
    a_prime = copy(a)
    x_control, x_target = a[control], a[target]
    z_control, z_target = a[n+control], a[n+target]
    a_prime[target] = (x_control + x_target) % 2
    a_prime[n+control] = (z_control + z_target) % 2
    r = ((x_control * z_target) * (x_target + z_control + 1)) % 2
    return a_prime, r
end

function cz_action_on_bsf(control, target, a)
    n = div(length(a), 2)
    a_prime = copy(a)
    x_control, x_target = a[control], a[target]
    z_control, z_target = a[n+control], a[n+target]
    a_prime[n+control] = (z_control + x_target) % 2
    a_prime[n+target] = (z_target + x_control) % 2
    r = ((x_control * x_target * z_control) + (x_control * x_target * z_target)) % 2
    return a_prime, r
end

function apply_pauli_basis_action(action_fun, q, hermitian, pauli_dict::PauliString)
    N = length(hermitian)
    hermitian_prime = zeros(Float64, N)
    for i in 1:N
        coeff = hermitian[i]
        pauli_vec = pauli_dict.int_to_bit[i]
        pauli_mapped, phase = action_fun(q, pauli_vec)
        index = pauli_dict.bit_to_int[pauli_mapped]
        hermitian_prime[index] = coeff * ((-1)^phase)
    end
    return hermitian_prime
end

function hadamard_action_on_pauli_basis(q, hermitian, pauli_dict::PauliString)
    return apply_pauli_basis_action(hadamard_action_on_bsf, q, hermitian, pauli_dict)
end

function phase_action_on_pauli_basis(q, hermitian, pauli_dict::PauliString)
    return apply_pauli_basis_action(phase_action_on_bsf, q, hermitian, pauli_dict)
end

function cx_action_on_pauli_basis(c, t, hermitian, pauli_dict::PauliString)
    action_fun(a) = cx_action_on_bsf(c, t, a)
    return apply_pauli_basis_action((q,a)->action_fun(a), nothing, hermitian, pauli_dict)
end

function cz_action_on_pauli_basis(c, t, hermitian, pauli_dict::PauliString)
    action_fun(a) = cz_action_on_bsf(c, t, a)
    return apply_pauli_basis_action((q,a)->action_fun(a), nothing, hermitian, pauli_dict)
end

# ------------------------------------------------------------------------------
# Direct construction of T^U (|+><+|)^⊗n T^{†U} in the Pauli basis
# ------------------------------------------------------------------------------

"""
Return the local Pauli coefficient for one qubit.

local_type ∈ (:I, :X, :Y, :Z)
is_red == true means apply T|+><+|T†
is_red == false means use |+><+|

Green:
    |+><+| has expectation coordinates (1, 1, 0, 0)

Red:
    T|+><+|T† has expectation coordinates (1, 1/√2, 1/√2, 0)
"""
function one_qubit_local_coeff(local_type::Symbol, is_red::Bool)
    if !is_red
        if local_type == :I || local_type == :X
            return 1.0
        else
            return 0.0
        end
    else
        if local_type == :I
            return 1.0
        elseif local_type == :X || local_type == :Y
            return 1 / sqrt(2)
        else
            return 0.0
        end
    end
end

"""
Decode the local Pauli type on qubit q from a BSF vector a.

Convention:
    (x,z) = (0,0) -> I
    (1,0) -> X
    (1,1) -> Y
    (0,1) -> Z
"""
function local_pauli_type(a, q, n)
    x = a[q]
    z = a[n + q]
    if x == 0 && z == 0
        return :I
    elseif x == 1 && z == 0
        return :X
    elseif x == 1 && z == 1
        return :Y
    else
        return :Z
    end
end

"""
Construct the Pauli vector of

    ρ_U = T^⊗U (|+><+|)^⊗n T^{†⊗U}

directly in the n-qubit Pauli basis.
"""
function initial_magic_product_state_pauli(n::Int, U::Vector{Int}, pauli_dict::PauliString)
    Uset = Set(U)
    coeffs = zeros(Float64, length(pauli_dict.bit_strings))

    for (i, a) in enumerate(pauli_dict.bit_strings)
        coeff = 1.0
        for q in 1:n
            local_type = local_pauli_type(a, q, n)
            coeff *= one_qubit_local_coeff(local_type, q in Uset)
            if coeff == 0.0
                break
            end
        end
        coeffs[i] = coeff
    end

    return coeffs
end

# ------------------------------------------------------------------------------
# Metadata helpers
# ------------------------------------------------------------------------------

function record_label(rec, k::Int)
    U = [Int(x) + 1 for x in rec["red"]]
    edges = [(Int(e[1]) + 1, Int(e[2]) + 1) for e in rec["edges"]]

    U_str = isempty(U) ? "" : join(U, ",")
    E_str = isempty(edges) ? "" : join(["($(i),$(j))" for (i,j) in edges], ",")

    return "state_$(k); U={$U_str}; E={$E_str}"
end

# ------------------------------------------------------------------------------
# Build the (G,U) resource projector in the Pauli basis from one JSON record
# ------------------------------------------------------------------------------

"""
Expected JSON record fields:
    n
    edges   # 0-based pairs from Python
    red     # 0-based red vertices from Python

Interpretation:
    red = U
"""
function resource_projector_pauli_from_record(rec, pauli_dict::PauliString)
    n = Int(rec["n"])
    U = [Int(x) + 1 for x in rec["red"]]

    rho = initial_magic_product_state_pauli(n, U, pauli_dict)

    for e in rec["edges"]
        i = Int(e[1]) + 1
        j = Int(e[2]) + 1
        rho = cz_action_on_pauli_basis(i, j, rho, pauli_dict)
    end

    return rho
end

# ------------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------------

function generate_resource_projectors_pauli(n::Int; save_hdf5::Bool=true)
    json_file = joinpath(data_folder, "bicolored_graphs_n$(n).json")
    records = JSON3.read(read(json_file, String))

    pauli_dict = PauliString(n)
    num_cases = length(records)
    dim_pauli = length(pauli_dict.bit_strings)

    println("Loaded $num_cases bicolored graphs from $json_file")

    M = Matrix{Float64}(undef, dim_pauli, num_cases)
    labels = Vector{String}(undef, num_cases)

    for (k, rec) in enumerate(records)
        M[:, k] = resource_projector_pauli_from_record(rec, pauli_dict)
        labels[k] = record_label(rec, k)
        println("Processed case $k / $num_cases")
    end

    if save_hdf5
        h5_file = joinpath(data_folder, "resource_projectors_pauli_n$(n).h5")
        h5open(h5_file, "w") do file
            file["data"] = M
            file["labels"] = labels
        end
        println("Saved Pauli-basis resource projectors to $h5_file")
    end

    return M, labels
end

# ------------------------------------------------------------------------------
# Run for fixed n
# ------------------------------------------------------------------------------

n = 6
generate_resource_projectors_pauli(n; save_hdf5=true)