################################################################################
# robustness_bicolored_resource_states.jl
#
# Load Pauli-basis vectors for the resource states
#
#     |G,U> = T^⊗n H^⊗U |G>
#
# and compute robustness values relative to a chosen phase-space matrix.
################################################################################

using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5

data_folder = joinpath(base_dir, "data")
local_robustness_folder = joinpath(base_dir, "local_robustness")
keys_folder = joinpath(base_dir, "keys")

push!(LOAD_PATH, local_robustness_folder)

include(joinpath(local_robustness_folder, "robustness_functions.jl"))

"""
Load an HDF5 dataset by name, falling back to the only dataset if needed.
"""
function load_hdf5_dataset(file_path::String; preferred_dataset::String="data")
    BicoloredRobustness.MaxWeightSeed.validate_maxweight_file(file_path)
    h5open(file_path, "r") do file
        ks = collect(keys(file))
        if preferred_dataset in ks
            return read(file[preferred_dataset])
        elseif length(ks) == 1
            return read(file[ks[1]])
        else
            error("Could not find dataset '$preferred_dataset' in $file_path. Available datasets: $(ks)")
        end
    end
end

"""
Load labels if present.
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

# phase space options include: "deterministic", "stabilizer_states", "cnc_operators", "max_weight_operators_p_2"
phase_space_name = "max_weight_operators_p_2"

for n in 2:3
    states_file = joinpath(data_folder, "bicolored_resource_projectors_pauli_n$(n).h5")
    phase_space_file = joinpath(data_folder, "$(phase_space_name)_$(n).h5")

    states = load_hdf5_dataset(states_file; preferred_dataset="data")
    labels = load_hdf5_labels(states_file)
    phase_space = load_hdf5_dataset(phase_space_file; preferred_dataset="data")

    println("Loaded states from: $states_file")
    println("Loaded phase space from: $phase_space_file")
    println("Number of states: ", size(states, 2))
    println("--------------------------------------------------")

    for j in 1:size(states, 2)
        state_label = labels === nothing ? "state_$(j)" : String(labels[j])
        println("State $j: $state_label")

        exp_name = "$(phase_space_name)_$(n)_bicolored_state_$(j)"
        robustness(
            phase_space,
            states[:, j];
            threshold=1e-16,
            save=true,
            file_loc=keys_folder,
            file_name=exp_name
        )
    end

    println("--------------------------------------------------")
end