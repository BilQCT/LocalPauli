################################################################################
# functions.jl
#
# This file defines common functions and utilities for generating graph states,
# applying quantum gate operations, converting states to the Pauli basis, and
# evaluating robustness. All functions include documentation comments.
################################################################################

using Pkg
# Activate the project environment located in "MyProject" folder
base_dir = dirname(@__DIR__)            # one level up from local_robustness
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5, JuMP, GLPK, LinearAlgebra, Combinatorics

# Set up directory paths relative to the project root
data_folder = joinpath(base_dir, "data")
libs_folder = joinpath(base_dir, "libs")
keys_folder = joinpath(base_dir, "keys")

# Optionally, add the libs folder to LOAD_PATH to allow includes:
push!(LOAD_PATH, libs_folder)

# Include the Pauli definitions from libs
include(joinpath(libs_folder, "pauli.jl"))

# ------------------------------------------------------------------------------
# Basic Gate Definitions
# ------------------------------------------------------------------------------

"""
    T, H, I2

Constants:
- `T`: 2×2 T gate matrix.
- `H`: 2×2 Hadamard matrix.
- `I2`: 2×2 Identity matrix.
"""
const T = Matrix{Complex}([1 0; 0 exp(im*π/4)])
const H = (1/sqrt(2)) * [1  1;
                         1 -1]
const I2 = Matrix{Complex}(I, 2, 2)

"""
    t_gate(n::Int, target::Int) -> Matrix{Complex}

Apply the single-qubit T gate on the `target` qubit of an n-qubit system.
Returns the full n-qubit gate as the Kronecker product of individual matrices.
"""
function t_gate(n::Int, target::Int)
    op_list = [I2 for _ in 1:n]
    op_list[target] = T
    return kron(op_list...)
end

"""
    hadamard_gate(n::Int, target::Int) -> Matrix{Complex}

Apply the single-qubit Hadamard gate on the `target` qubit of an n-qubit system.
Returns the full n-qubit gate as the Kronecker product.
"""
function hadamard_gate(n::Int, target::Int)
    op_list = [I2 for _ in 1:n]
    op_list[target] = H
    return kron(op_list...)
end

"""
    plus_state(n::Int) -> Vector{Complex}

Construct the tensor product of n copies of the |+> state (with |+> = [1, 1]^T).
"""
function plus_state(n::Int)
    plus = Vector{Complex}([1; 1])
    plus_n = [1]
    for _ in 1:n
        plus_n = kron(plus_n, plus)
    end
    return plus_n
end

"""
    controlled_z_gate(n::Int, control::Int, target::Int) -> Matrix{Float64}

Generate the Controlled-Z gate for an n-qubit system acting on the specified
`control` and `target` qubits. The gate applies a phase flip (i.e. -1) when both
qubits are in the |1⟩ state.
"""
function controlled_z_gate(n::Int, control::Int, target::Int)
    dim = 2^n
    CZ = Matrix{Float64}(I, dim, dim)
    for i in 0:dim-1
        # Get binary representation with n digits (left-to-right order)
        binary_repr = reverse(digits(i, base=2, pad=n))
        if binary_repr[control] == 1 && binary_repr[target] == 1
            CZ[i+1, i+1] = -1
        end
    end
    return CZ
end

# ------------------------------------------------------------------------------
# Graph State Generation (State-Vector Approach)
# ------------------------------------------------------------------------------

"""
    generate_graph_state_vectors(n::Int) -> (Vector{Vector{Complex}}, Vector{String})

Generates three types of n-qubit graph states (line, cycle, complete) using the
|+⟩^n state and applying Controlled-Z gates appropriately. Returns a tuple:
- A vector of state vectors.
- A vector of corresponding names.
"""
function generate_graph_state_vectors(n::Int)
    Plus = plus_state(n)
    
    # Line graph: CZ between adjacent qubits
    Psi_Ln = copy(Plus)
    for i in 1:n-1
        Psi_Ln = controlled_z_gate(n, i, i+1) * Psi_Ln
    end
    
    # Cycle graph: line graph with an extra CZ between qubit 1 and n
    Psi_Cn = copy(Plus)
    for i in 1:n
        if i < n
            Psi_Cn = controlled_z_gate(n, i, i+1) * Psi_Cn
        else
            Psi_Cn = controlled_z_gate(n, 1, n) * Psi_Cn
        end
    end
    
    # Complete graph: apply CZ for every pair of qubits
    pairs = [p for p in combinations(1:n, 2)]
    Psi_Kn = copy(Plus)
    for p in pairs
        Psi_Kn = controlled_z_gate(n, p[1], p[2]) * Psi_Kn
    end
    
    return [Psi_Ln, Psi_Cn, Psi_Kn], ["L$n", "C$n", "K$n"]
end

# ------------------------------------------------------------------------------
# Pauli Basis and Tableau Actions
# ------------------------------------------------------------------------------

# The following assumes that you have a module or definitions for the Pauli
# basis (e.g., a type `PauliString` and function `pauli`). Here we include the
# external file that defines these.
include("../libs/pauli.jl")

"""
    convert_matrix_to_pauli_basis(n::Int, A) -> Vector{Float64}

Converts a density matrix (or operator) A into its expansion coefficients in
the n-qubit Pauli basis.
"""
function convert_matrix_to_pauli_basis(n::Int, A)
    pauli_dict = PauliString(n)
    vector_A = Float64[]
    for a in pauli_dict.bit_strings
        Ta = pauli(a)
        A_a = real(tr(Ta * A) / (2^n))
        push!(vector_A, A_a)
    end
    return vector_A
end

# --- Tableau Action Functions ---

"""
    hadamard_action_on_bsf(q::Int, a::Vector{Int}) -> (Vector{Int}, Int)

Simulates the effect of a Hadamard gate on the binary symplectic (tableau)
representation of a Pauli operator. Returns the transformed tableau and a phase.
"""
function hadamard_action_on_bsf(q::Int, a::Vector{Int})
    n = div(length(a), 2)
    a_prime = copy(a)
    a_prime[q] = a[n+q]
    a_prime[n+q] = a[q]
    r = (a[q] * a[n+q]) % 2
    return a_prime, r
end

"""
    phase_action_on_bsf(q::Int, a::Vector{Int}) -> (Vector{Int}, Int)

Simulates the effect of a phase gate on the tableau representation of a Pauli
operator. Returns the transformed tableau and a phase.
"""
function phase_action_on_bsf(q::Int, a::Vector{Int})
    n = div(length(a), 2)
    a_prime = copy(a)
    a_prime[n+q] = (a[q] + a[n+q]) % 2
    r = (a[q] * a[n+q]) % 2
    return a_prime, r
end

"""
    cx_action_on_bsf(control::Int, target::Int, a::Vector{Int}) -> (Vector{Int}, Int)

Simulates the action of a CNOT gate (control-target) on the tableau representation.
Returns the transformed tableau and a phase.
"""
function cx_action_on_bsf(control::Int, target::Int, a::Vector{Int})
    n = div(length(a), 2)
    a_prime = copy(a)
    x_control, x_target = a[control], a[target]
    z_control, z_target = a[n+control], a[n+target]
    a_prime[target] = (x_control + x_target) % 2
    a_prime[n+control] = (z_control + z_target) % 2
    r = ((x_control * z_target) * (x_target + z_control + 1)) % 2
    return a_prime, r
end

"""
    cz_action_on_bsf(control::Int, target::Int, a::Vector{Int}) -> (Vector{Int}, Int)

Simulates the action of a Controlled-Z gate on the tableau representation.
Returns the transformed tableau and a phase.
"""
function cz_action_on_bsf(control::Int, target::Int, a::Vector{Int})
    n = div(length(a), 2)
    a_prime = copy(a)
    x_control, x_target = a[control], a[target]
    z_control, z_target = a[n+control], a[n+target]
    a_prime[n+control] = (z_control + x_target) % 2
    a_prime[n+target] = (z_target + x_control) % 2
    r = ((x_control * x_target * z_control) + (x_control * x_target * z_target)) % 2
    return a_prime, r
end

"""
    apply_pauli_basis_action(action_fun, q, hermitian::Vector{Real}, pauli_dict::PauliString)
       -> Vector{Real}

Helper function that applies a given tableau action (e.g., Hadamard, Phase)
to each element of the Pauli basis expansion of a Hermitian operator.
"""
function apply_pauli_basis_action(action_fun, q, hermitian::Vector{Real}, pauli_dict::PauliString)
    N = length(hermitian)
    hermitian_prime = ones(Float64, N)
    for i in 1:N
        coeff = hermitian[i]
        pauli_vec = pauli_dict.int_to_bit[i]
        pauli_mapped, phase = action_fun(q, pauli_vec)
        index = pauli_dict.bit_to_int[pauli_mapped]
        hermitian_prime[index] = coeff * ((-1)^phase)
    end
    return hermitian_prime
end

# Specific wrappers for common operations on the Pauli basis:
function hadamard_action_on_pauli_basis(q, hermitian::Vector{Real}, pauli_dict::PauliString)
    return apply_pauli_basis_action(hadamard_action_on_bsf, q, hermitian, pauli_dict)
end

function phase_action_on_pauli_basis(q, hermitian::Vector{Real}, pauli_dict::PauliString)
    return apply_pauli_basis_action(phase_action_on_bsf, q, hermitian, pauli_dict)
end

function cx_action_on_pauli_basis(c, t, hermitian::Vector{Real}, pauli_dict::PauliString)
    action_fun(a) = cx_action_on_bsf(c, t, a)
    return apply_pauli_basis_action((q,a)->action_fun(a), nothing, hermitian, pauli_dict)
end

function cz_action_on_pauli_basis(c, t, hermitian::Vector{Real}, pauli_dict::PauliString)
    action_fun(a) = cz_action_on_bsf(c, t, a)
    return apply_pauli_basis_action((q,a)->action_fun(a), nothing, hermitian, pauli_dict)
end

# ------------------------------------------------------------------------------
# Graph State Generation (Pauli Basis Approach)
# ------------------------------------------------------------------------------

"""
    generate_graph_states(magic_state, n, pauli_dict::PauliString)
        -> (Vector{Vector{Real}}, Vector{String})

Generates graph states starting from a magic state (in the Pauli basis)
for an n-qubit system. It returns three types of graph states (line, cycle,
complete) along with their names.
"""
function generate_graph_states(magic_state, n, pauli_dict::PauliString)
    # Line graph
    Ln = copy(magic_state)
    for i in 1:n-1
        Ln = cz_action_on_pauli_basis(i, i+1, Ln, pauli_dict)
    end
    # Cycle graph
    Cn = copy(magic_state)
    for i in 1:n
        if i < n
            Cn = cz_action_on_pauli_basis(i, i+1, Cn, pauli_dict)
        else
            Cn = cz_action_on_pauli_basis(1, n, Cn, pauli_dict)
        end
    end
    # Complete graph
    pairs = [p for p in combinations(1:n, 2)]
    Kn = copy(magic_state)
    for p in pairs
        Kn = cz_action_on_pauli_basis(p[1], p[2], Kn, pauli_dict)
    end
    return [Ln, Cn, Kn], ["L$n", "C$n", "K$n"]
end


# ------------------------------------------------------------------------------
# Robustness Evaluation
# ------------------------------------------------------------------------------

function robustness(input_array, target_state; threshold=1e-16, save=true, file_loc::String=keys_folder, file_name="result")
    model = Model(GLPK.Optimizer)
    N = size(input_array, 2)
    @variable(model, x[1:N])
    @variable(model, u[1:N] ≥ 0)
    @objective(model, Min, sum(u))
    @constraint(model, [i=1:N], u[i] ≥ x[i])
    @constraint(model, [i=1:N], u[i] ≥ -x[i])
    @constraint(model, input_array * x .== target_state)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL
        error("Optimization failed with status: ", termination_status(model))
    end
    x_value = value.(x)
    obj_value = objective_value(model)
    println("Optimal objective value (min ||x||₁): ", obj_value)
    significant_terms = findall(x -> abs(x) > threshold, x_value)
    decomposition_array = vcat(transpose(x_value[significant_terms]), input_array[:, significant_terms])
    if save
        h5_file = joinpath(file_loc, file_name * ".h5")
        h5open(h5_file, "w") do file
            file[file_name] = decomposition_array
        end
    end
end

# ------------------------------------------------------------------------------
# Additional Utility Functions
# ------------------------------------------------------------------------------

"""
    t_state(n::Int) -> Vector{Real}

Generates a magic state for n qubits based on a predefined pattern.
"""
function t_state(n::Int)
    t = [1, 1/sqrt(2), 1/sqrt(2), 0]
    T_state = [1]
    for _ in 1:n
        T_state = kron(T_state, t)
    end
    return Vector{Real}(T_state)
end

"""
    all_dit_strings(base, n) -> Vector{Vector{Int}}

Generates all base-`base` binary strings of length `n`. (This is a dummy
implementation; replace with your actual function if needed.)
"""
function all_dit_strings(base, n)
    return [digits(i, base=base, pad=n) for i in 0:(base^n - 1)]
end
