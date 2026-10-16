# Historical E(G) T^U |+> resource convention; NOT the manuscript robustness target.
# Require an explicit opt-in to prevent accidentally generating publication data.
if abspath(PROGRAM_FILE) == @__FILE__ && !("--legacy" in ARGS)
    error("Historical state convention. Use the bicolored manuscript pipeline, or pass --legacy explicitly.")
end
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

for n in 2:4
    states_file = joinpath(data_folder, "resource_projectors_pauli_n$(n).h5")
    phase_space_file = joinpath(data_folder, "deterministic_$(n).h5")

    global states = nothing
    global labels = nothing
    global phase_space = nothing

    h5open(states_file, "r") do file
        global states = read(file["data"])
        if "labels" in keys(file)
            global labels = read(file["labels"])
        end
    end

    h5open(phase_space_file, "r") do file
        global phase_space = read(file["data"])
    end

    println("Loaded states from: $states_file")
    println("Loaded phase space from: $phase_space_file")
    println("Number of states: ", size(states, 2))
    println("--------------------------------------------------")

    for j in 1:size(states, 2)
        if labels === nothing
            println("State $j:")
        else
            println("State $j: $(labels[j])")
        end

        exp_name = "deterministic" * "_" * "n$(n)" * "_" * "state_$(j)"
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