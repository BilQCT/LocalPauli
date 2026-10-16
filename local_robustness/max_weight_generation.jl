using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5
using LinearAlgebra
using Combinatorics

# Load local Clifford / Bloch-vector machinery.
include(joinpath(@__DIR__, "..", "libs", "LocalCliffordBloch.jl"))
using .LocalCliffordBloch

# -----------------------------------------------------------------------------
# Efficient max-weight / NBL fiducials
# Map 0,1,2 in Z_3 to x,y,z on each qubit.
# -----------------------------------------------------------------------------

# Seed construction is shared with the dictionary loader and regression tests.
# See libs/MaxWeightSeed.jl for the explicit modular arithmetic.

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
                attributes(file)["seed_convention"] = LocalCliffordBloch.MaxWeightSeed.SEED_CONVENTION
                attributes(file)["orbit_convention"] = LocalCliffordBloch.MaxWeightSeed.ORBIT_CONVENTION
                attributes(file)["pauli_convention"] = "lexicographic I,X,Y,Z; Tr(P A)"
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

if abspath(PROGRAM_FILE) == @__FILE__
    ns = isempty(ARGS) ? (2:4) : parse.(Int, ARGS)
    save_max_weight_orbits_by_p_h5(ns; data_folder=data_folder)
end
