################################################################################
# bicolored_resource_state_generation.jl
#
# Generate n-qubit resource states
#
#     |G,U> = T^⊗n H^⊗U |G>
#
# in the unnormalized Pauli/Bloch convention used in your previous scripts.
#
# Output columns are Pauli vectors of the projectors |G,U><G,U|.
################################################################################

using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using JSON3, HDF5, LinearAlgebra

data_folder = joinpath(base_dir, "data")
libs_folder = joinpath(base_dir, "libs")

push!(LOAD_PATH, libs_folder)
include(joinpath(libs_folder, "pauli.jl"))


function all_dit_strings(base::Int, n::Int)
    return [digits(i, base=base, pad=n) for i in 0:(base^n - 1)]
end

function print_centered_message(msg::String, width::Int=80, fill_char::String="-")
    padding = max(0, div(width - length(msg), 2))
    println(rpad(lpad(msg, padding + length(msg), fill_char), width, fill_char))
end

# ------------------------------------------------------------------------------
# Clifford actions on BSF
# ------------------------------------------------------------------------------

function hadamard_action_on_bsf(q::Int, a::Vector{Int})
    n = div(length(a), 2)
    a_prime = copy(a)
    a_prime[q] = a[n + q]
    a_prime[n + q] = a[q]
    r = (a[q] * a[n + q]) % 2
    return a_prime, r
end

function cz_action_on_bsf(control::Int, target::Int, a::Vector{Int})
    n = div(length(a), 2)
    a_prime = copy(a)

    x_control, x_target = a[control], a[target]
    z_control, z_target = a[n + control], a[n + target]

    a_prime[n + control] = (z_control + x_target) % 2
    a_prime[n + target]  = (z_target + x_control) % 2

    r = (x_control * x_target * z_control + x_control * x_target * z_target) % 2
    return a_prime, r
end

function apply_pauli_basis_action(action_fun, q, hermitian::Vector{Float64}, pauli_dict::PauliString)
    N = length(hermitian)
    hermitian_prime = zeros(Float64, N)

    for i in 1:N
        coeff = hermitian[i]
        coeff == 0.0 && continue

        pauli_vec = pauli_dict.int_to_bit[i]
        pauli_mapped, phase = action_fun(q, pauli_vec)
        index = pauli_dict.bit_to_int[pauli_mapped]
        hermitian_prime[index] += coeff * ((-1)^phase)
    end

    return hermitian_prime
end

function hadamard_action_on_pauli_basis(q::Int, hermitian::Vector{Float64}, pauli_dict::PauliString)
    return apply_pauli_basis_action(hadamard_action_on_bsf, q, hermitian, pauli_dict)
end

function cz_action_on_pauli_basis(c::Int, t::Int, hermitian::Vector{Float64}, pauli_dict::PauliString)
    action_fun(a) = cz_action_on_bsf(c, t, a)
    return apply_pauli_basis_action((q, a) -> action_fun(a), nothing, hermitian, pauli_dict)
end

# ------------------------------------------------------------------------------
# T action on the Pauli basis
# ------------------------------------------------------------------------------

function local_pauli_type(a::Vector{Int}, q::Int, n::Int)
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

function set_local_pauli_type(a::Vector{Int}, q::Int, n::Int, typ::Symbol)
    b = copy(a)

    if typ == :I
        b[q] = 0
        b[n + q] = 0
    elseif typ == :X
        b[q] = 1
        b[n + q] = 0
    elseif typ == :Y
        b[q] = 1
        b[n + q] = 1
    elseif typ == :Z
        b[q] = 0
        b[n + q] = 1
    else
        error("Unknown Pauli type: $typ")
    end

    return b
end

function t_action_on_pauli_basis(q::Int, coeffs::Vector{Float64}, pauli_dict::PauliString)
    n = div(length(pauli_dict.bit_strings[1]), 2)
    out = zeros(Float64, length(coeffs))
    invsqrt2 = 1 / sqrt(2)

    for i in eachindex(coeffs)
        c = coeffs[i]
        c == 0.0 && continue

        a = pauli_dict.int_to_bit[i]
        typ = local_pauli_type(a, q, n)

        if typ == :I || typ == :Z
            out[i] += c

        elseif typ == :X
            # X -> (X + Y)/sqrt(2)
            out[i] += c * invsqrt2
            aY = set_local_pauli_type(a, q, n, :Y)
            idxY = pauli_dict.bit_to_int[aY]
            out[idxY] += c * invsqrt2

        elseif typ == :Y
            # Y -> (Y - X)/sqrt(2)
            out[i] += c * invsqrt2
            aX = set_local_pauli_type(a, q, n, :X)
            idxX = pauli_dict.bit_to_int[aX]
            out[idxX] -= c * invsqrt2
        end
    end

    return out
end

function global_t_action_on_pauli_basis(coeffs::Vector{Float64}, pauli_dict::PauliString)
    n = div(length(pauli_dict.bit_strings[1]), 2)
    out = copy(coeffs)
    for q in 1:n
        out = t_action_on_pauli_basis(q, out, pauli_dict)
    end
    return out
end

# ------------------------------------------------------------------------------
# Initial state in the user's unnormalized convention
# ------------------------------------------------------------------------------

"""
Unnormalized Pauli vector for |+><+|^⊗n.

For one qubit:
    |+><+|  ↦  (1, 1, 0, 0)

in the user's convention.
"""
function plus_state_pauli(n::Int)
    plus_pauli = [1.0, 1.0, 0.0, 0.0]
    out = [1.0]
    for _ in 1:n
        out = kron(out, plus_pauli)
    end
    return Vector{Float64}(out)
end

function resource_projector_pauli(n::Int,
                                  edges::Vector{Tuple{Int,Int}},
                                  U::Vector{Int},
                                  pauli_dict::PauliString)
    rho = plus_state_pauli(n)

    for (i, j) in edges
        rho = cz_action_on_pauli_basis(i, j, rho, pauli_dict)
    end

    for q in U
        rho = hadamard_action_on_pauli_basis(q, rho, pauli_dict)
    end

    rho = global_t_action_on_pauli_basis(rho, pauli_dict)

    return rho
end

# ------------------------------------------------------------------------------
# JSON handling
# ------------------------------------------------------------------------------

function load_bicolored_graph_records(n::Int)
    json_file = joinpath(data_folder, "bicolored_graphs_n$(n).json")
    records = JSON3.read(read(json_file, String))
    return records, json_file
end

function parse_record(rec)
    n = Int(rec["n"])
    edges = Tuple{Int,Int}[(Int(e[1]) + 1, Int(e[2]) + 1) for e in rec["edges"]]
    U = [Int(x) + 1 for x in rec["red"]]
    name = haskey(rec, "name") ? String(rec["name"]) : ""
    return n, edges, U, name
end

# ------------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------------

function generate_resource_projectors_pauli(n::Int; save_hdf5::Bool=true)
    records, json_file = load_bicolored_graph_records(n)
    pauli_dict = PauliString(n)

    num_cases = length(records)
    dim_pauli = length(pauli_dict.bit_strings)

    println("\n" * "="^80)
    print_centered_message(" GENERATING BICOLORED RESOURCE PROJECTORS, n = $n ", 80, "=")
    println("="^80)
    println("Loaded $num_cases graph/color records from:")
    println("  $json_file")
    println()

    M = Matrix{Float64}(undef, dim_pauli, num_cases)
    labels = Vector{String}(undef, num_cases)

    for (k, rec) in enumerate(records)
        n_rec, edges, U, name = parse_record(rec)

        if n_rec != n
            error("Record $k has n = $n_rec but expected n = $n.")
        end

        M[:, k] = resource_projector_pauli(n, edges, U, pauli_dict)

        edge_str = join(["($(i),$(j))" for (i, j) in edges], ",")
        U_str = join(U, ",")
        labels[k] = isempty(name) ? "state_$(k); U={$U_str}; E={$edge_str}" : name

        println("Processed state $k / $num_cases")
    end

    if save_hdf5
        h5_file = joinpath(data_folder, "bicolored_resource_projectors_pauli_n$(n).h5")
        h5open(h5_file, "w") do file
            attributes(file)["state_family"] = "T^n H^U E(G) |+>^n"
            attributes(file)["pauli_convention"] = "lexicographic I,X,Y,Z; Tr(P rho)"
            file["data"] = M
            file["labels"] = labels
        end
        println()
        println("Saved Pauli-basis resource projectors to:")
        println("  $h5_file")
    end

    return M, labels
end

if abspath(PROGRAM_FILE) == @__FILE__
    for n in (isempty(ARGS) ? (2:4) : parse.(Int, ARGS))
        generate_resource_projectors_pauli(n; save_hdf5=true)
    end
end