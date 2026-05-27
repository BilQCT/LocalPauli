using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5
using LinearAlgebra
using Combinatorics

# Load local Clifford / Bloch-vector machinery.
include("../libs/LocalCliffordBloch.jl")
using .LocalCliffordBloch

# -----------------------------------------------------------------------------
# Efficient max-weight / NBL fiducials
# Map 0,1,2 in Z_3 to x,y,z on each qubit.
# -----------------------------------------------------------------------------

function z3_to_local_symplectic(u::Vector{Int})
    n = length(u)
    a = zeros(Int, 2*n)

    for i in 1:n
        if u[i] == 0       # X_i
            a[i] = 1
        elseif u[i] == 1   # Y_i = X_i Z_i
            a[i] = 1
            a[n+i] = 1
        elseif u[i] == 2   # Z_i
            a[n+i] = 1
        else
            error("Entries of u must lie in {0,1,2}.")
        end
    end

    return a
end

function gamma_elementary_p(u::Vector{Int}, p::Int)
    n = length(u)

    if p < 1
        error("p must be positive.")
    elseif p > n
        return 0
    end

    total_z3 = 0
    for I in combinations(1:n, p)
        prod_z3 = 1
        for i in I
            prod_z3 = mod(prod_z3 * u[i], 3)
        end
        total_z3 = mod(total_z3 + prod_z3, 3)
    end

    return mod(total_z3, 2)
end

function max_weight_bloch_vector(n::Int, p::Int; pauli_info::PauliString=PauliString(n))
    if n < 2
        error("Need n ≥ 2.")
    elseif p <= 1
        error("Need p > 1.")
    elseif p > n
        error("Need p ≤ n for this elementary degree-p fiducial.")
    end

    C = zeros(Int, 4^n)
    bit_to_int = pauli_info.bit_to_int

    C[bit_to_int[zeros(Int, 2*n)]] = 1

    for u_tuple in Iterators.product(ntuple(_ -> 0:2, n)...)
        u = collect(u_tuple)
        a = z3_to_local_symplectic(u)
        γ = gamma_elementary_p(u, p)
        C[bit_to_int[a]] = (-1)^γ
    end

    return C
end

function all_max_weight_fiducials(n::Int; pauli_info::PauliString=PauliString(n))
    return Dict(p => max_weight_bloch_vector(n, p; pauli_info=pauli_info) for p in 2:n)
end


# -----------------------------------------------------------------------------
# Orbit generation and HDF5 saving
# -----------------------------------------------------------------------------
#
# Output convention:
#   data[:, j] is the j-th max-weight operator in Bloch/Pauli coordinates.
#
# The file also stores:
#   labels      :: string labels, e.g. "n=4,p=2,orbit=17"
#   p_labels    :: integer p value for each column
#   orbit_index :: integer index within the p-orbit for each column
#
# The important key for compatibility with your existing loader is:
#   "data"
# -----------------------------------------------------------------------------

function orbit_matrix_from_fiducial(CG::LocalCliffordGroup, C::AbstractVector)
    orbit = local_clifford_orbit_of_point(CG, collect(C))

    if isempty(orbit)
        error("The local Clifford orbit is empty. Something is wrong.")
    end

    # columns are operators
    M = hcat(orbit...)

    # The max-weight coefficients should be in {-1,0,1}, so Int8 is enough.
    return Int8.(M)
end

function save_max_weight_orbits_by_p_h5(ns; data_folder="./data")
    mkpath(data_folder)

    for n in ns
        println("Generating LocalCliffordGroup($n)...")
        CG = LocalCliffordGroup(n)
        pauli_info = CG.PauliInfo

        for p in 2:2
            println("  Generating max-weight orbit for n = $n, p = $p")

            C = max_weight_bloch_vector(n, p; pauli_info=pauli_info)
            orbit = local_clifford_orbit_of_point(CG, C)

            data = Int8.(hcat(orbit...))

            filename = joinpath(data_folder, "max_weight_operators_p_$(p)_$(n).h5")

            h5open(filename, "w") do file
                file["data"] = data
                file["pauli_labels"] = pauli_info.pauli_strings
                file["p"] = p
                file["n"] = n
                file["fiducial"] = Int8.(C)
            end

            println("    saved: $filename")
            println("    size(data): ", size(data))
        end
    end
end


# -----------------------------------------------------------------------------
# Run as a script
# -----------------------------------------------------------------------------
#
# Edit this range as needed. For example:
#   ns_to_generate = 2:3
#   ns_to_generate = 2:4
# -----------------------------------------------------------------------------

data_folder = joinpath(base_dir, "data")

save_max_weight_orbits_by_p_h5(4:4; data_folder=data_folder)
