using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5
using LinearAlgebra

# Load your local Clifford / Bloch-vector machinery.
# Adjust this path as needed.
include("LocalCliffordBloch.jl")

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

function kron_power(v::AbstractVector, n::Int)
    if n < 1
        error("n must be positive.")
    end

    out = copy(v)
    for _ in 2:n
        out = kron(out, v)
    end

    return out
end


function save_hdf5_dataset(filename::AbstractString, data; labels=nothing)
    mkpath(dirname(filename))

    h5open(filename, "w") do file
        file["data"] = data

        if labels !== nothing
            file["labels"] = labels
        end
    end
end


"""
    generate_local_stabilizer_states(n; data_folder="./data")

Generate local stabilizer states in Bloch-vector form as the local Clifford orbit
of |+><+|^{⊗ n}.

Single-qubit |+><+| has Bloch vector

    [1, 1, 0, 0]

in the convention

    [I, X, Y, Z].

The output matrix has size

    4^n × 6^n

with one local stabilizer state per column.
"""

using .LocalCliffordBloch

function generate_local_stabilizer_states(n::Int; data_folder="./data")
    if n < 1
        error("n must be positive.")
    end

    CG = LocalCliffordGroup(n)

    # Single-qubit |+><+| = (I + X)/2.
    P1 = Vector{Int}([1, 1, 0, 0])

    # Fiducial product local stabilizer state.
    Pn = kron_power(P1, n)

    # Orbit under local Clifford group.
    LS_orbit = local_clifford_orbit_of_point(CG, Pn)

    # Store as matrix: columns are operators/states.
    LS = hcat(LS_orbit...)

    # Use a compact integer type if coefficients are only -1,0,1.
    LS = Int8.(LS)

    # Optional Pauli labels for rows.
    labels = CG.PauliInfo.pauli_strings

    outfile = joinpath(data_folder, "local_stabilizer_states_$(n).h5")

    save_hdf5_dataset(outfile, LS; labels=labels)

    println("Saved local stabilizer states for n = $n")
    println("  file: ", outfile)
    println("  size(data): ", size(LS))
    println("  number of states/operators: ", size(LS, 2))
    println("  expected number: ", 6^n)

    return LS
end


# ----------------------------------------------------------------------
# Generate for desired n
# ----------------------------------------------------------------------

data_folder = "./data"

for n in 2:4
    generate_local_stabilizer_states(n; data_folder=data_folder)
end