################################################################################
# robustness_bicolored_resource_states_combined_phase_spaces.jl
#
# Load Pauli-basis vectors for the resource states
#
#     |G,U> = T^⊗n H^⊗U |G>
#
# and compute robustness values relative to either:
#
#   (i)  a single phase-space matrix, e.g.
#        "deterministic"
#
#   (ii) a combined phase-space matrix, e.g.
#        ["deterministic", "max_weight_operators_p_2"]
#
# Each phase-space file is expected to live in `data_folder` and have name
#
#     $(phase_space_name)_$(n).h5
#
# with HDF5 dataset key
#
#     "data"
#
# The combined phase space is formed by horizontal concatenation of columns.
################################################################################

using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5
using LinearAlgebra

data_folder = joinpath(base_dir, "data")
local_robustness_folder = joinpath(base_dir, "local_robustness")
keys_folder = joinpath(base_dir, "keys")

push!(LOAD_PATH, local_robustness_folder)

include(joinpath(local_robustness_folder, "robustness_functions.jl"))

# ------------------------------------------------------------------------------
# HDF5 loading helpers
# ------------------------------------------------------------------------------

"""
    load_hdf5_dataset(file_path; preferred_dataset="data")

Load an HDF5 dataset by name, falling back to the only dataset if needed.
"""
function load_hdf5_dataset(file_path::String; preferred_dataset::String="data")
    h5open(file_path, "r") do file
        ks = collect(keys(file))

        if preferred_dataset in ks
            return read(file[preferred_dataset])
        elseif length(ks) == 1
            return read(file[ks[1]])
        else
            error(
                "Could not find dataset '$preferred_dataset' in $file_path. " *
                "Available datasets: $(ks)"
            )
        end
    end
end


"""
    load_hdf5_labels(file_path)

Load labels if present. Returns `nothing` otherwise.
"""
function load_hdf5_labels(file_path::String)
    h5open(file_path, "r") do file
        if "labels" in keys(file)
            return read(file["labels"])
        else
            return nothing
        end
    end
end

# ------------------------------------------------------------------------------
# Phase-space combination helpers
# ------------------------------------------------------------------------------

"""
    phase_space_label(phase_space_spec)

Return a canonical string label for a single or combined phase space.

Examples:
    phase_space_label("deterministic")
        == "deterministic"

    phase_space_label(["deterministic", "max_weight_operators_p_2"])
        == "deterministic_plus_max_weight_operators_p_2"
"""
function phase_space_label(phase_space_spec)
    if phase_space_spec isa AbstractString
        return String(phase_space_spec)
    else
        return join(string.(phase_space_spec), "_plus_")
    end
end


"""
    load_single_phase_space(data_folder, phase_space_name, n)

Load one phase-space matrix from

    (data_folder)/(phase_space_name)_(n).h5

using dataset key "data".
"""
function load_single_phase_space(
    data_folder::String,
    phase_space_name::AbstractString,
    n::Int,
)
    phase_space_file = joinpath(data_folder, "$(phase_space_name)_$(n).h5")

    if !isfile(phase_space_file)
        error("Phase-space file does not exist: $phase_space_file")
    end

    phase_space = load_hdf5_dataset(phase_space_file; preferred_dataset="data")

    println("Loaded phase space from: $phase_space_file")
    println("  size: ", size(phase_space))

    return phase_space
end


"""
    unique_phase_space_columns(A)

Remove duplicate columns from a phase-space matrix.

This is useful when combining phase spaces, e.g. deterministic plus
max-weight, or max-weight orbits for several p-values.
"""
function unique_phase_space_columns(A::AbstractMatrix)
    seen = Set{Vector{eltype(A)}}()
    keep_cols = Int[]

    for j in axes(A, 2)
        col = Vector(A[:, j])

        if !(col in seen)
            push!(seen, col)
            push!(keep_cols, j)
        end
    end

    return A[:, keep_cols]
end


"""
    load_phase_space(data_folder, phase_space_spec, n; deduplicate=true)

Load a phase-space matrix.

`phase_space_spec` can be either:

    "deterministic"

or a vector such as:

    ["deterministic", "max_weight_operators_p_2"]

In the vector case, the corresponding phase-space matrices are horizontally
concatenated.
"""
function load_phase_space(
    data_folder::String,
    phase_space_spec,
    n::Int;
    deduplicate::Bool=true,
)
    if phase_space_spec isa AbstractString
        return load_single_phase_space(data_folder, phase_space_spec, n)
    end

    blocks = Matrix[]

    for name in phase_space_spec
        block = load_single_phase_space(data_folder, name, n)
        push!(blocks, block)
    end

    # Check row compatibility before hcat.
    row_counts = [size(block, 1) for block in blocks]
    if length(unique(row_counts)) != 1
        error(
            "Cannot combine phase spaces with different row counts: " *
            "$(collect(zip(phase_space_spec, row_counts)))"
        )
    end

    phase_space = hcat(blocks...)

    if deduplicate
        old_cols = size(phase_space, 2)
        phase_space = unique_phase_space_columns(phase_space)
        new_cols = size(phase_space, 2)
        println("Deduplicated columns: $old_cols -> $new_cols")
    end

    println("Combined phase space: ", phase_space_label(phase_space_spec))
    println("  combined size: ", size(phase_space))

    return phase_space
end

# ------------------------------------------------------------------------------
# Robustness run configuration
# ------------------------------------------------------------------------------

# Single phase-space examples:
# phase_space_spec = "deterministic"
# phase_space_spec = "cnc_operators"
# phase_space_spec = "max_weight_operators_p_2"

# Default combined phase space:
phase_space_spec = ["deterministic", "cnc_operators"]

# Later examples:
# phase_space_spec = ["deterministic", "max_weight_operators_p_2", "max_weight_operators_p_3"]
# phase_space_spec = ["deterministic", "local_stabilizer_states", "max_weight_operators_p_2"]

phase_space_name = phase_space_label(phase_space_spec)

# Choose the n-values to process.
n_values = 3:3

# Whether to remove duplicate columns after combining phase spaces.
deduplicate_phase_space = true

# Rank calculation can be useful, but can be expensive for larger matrices.
compute_phase_space_rank = true

# Robustness LP threshold.
robustness_threshold = 1e-16

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

for n in n_values
    println("==================================================")
    println("n = $n")
    println("Phase-space spec: ", phase_space_spec)
    println("Phase-space label: ", phase_space_name)
    println("==================================================")

    states_file = joinpath(data_folder, "bicolored_resource_projectors_pauli_n$(n).h5")

    if !isfile(states_file)
        error("Resource-state file does not exist: $states_file")
    end

    states = load_hdf5_dataset(states_file; preferred_dataset="data")
    labels = load_hdf5_labels(states_file)

    phase_space = load_phase_space(
        data_folder,
        phase_space_spec,
        n;
        deduplicate=deduplicate_phase_space,
    )

    println("Loaded states from: $states_file")
    println("Number of states: ", size(states, 2))
    println("State vector dimension: ", size(states, 1))
    println("Phase-space size: ", size(phase_space))

    if size(states, 1) != size(phase_space, 1)
        error(
            "Dimension mismatch for n = $n: " *
            "states have $(size(states, 1)) rows, " *
            "phase space has $(size(phase_space, 1)) rows."
        )
    end

    if compute_phase_space_rank
        println("Rank of phase space: ", rank(Matrix{Float64}(phase_space)))
    end

    println("--------------------------------------------------")

    for j in 1:size(states, 2)
        state_label = labels === nothing ? "state_$(j)" : String(labels[j])
        println("State $j: $state_label")

        exp_name = "$(phase_space_name)_$(n)_bicolored_state_$(j)"

        robustness(
            phase_space,
            states[:, j];
            threshold=robustness_threshold,
            save=true,
            file_loc=keys_folder,
            file_name=exp_name,
        )
    end

    println("--------------------------------------------------")
end
