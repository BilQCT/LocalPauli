################################################################################
# run_experiment.jl
#
# This driver script runs experiments by loading an input phase-space matrix,
# generating graph states, applying Hadamard and T operations in different orders,
# and evaluating the robustness of the resulting states.
#
# Two modes are supported:
#
#   :mbqc       -> These are resource states arising from a reduction technique used in mbqc.
#                       Graph states are generated first (from |+⟩^n). Then, 
#                       Hadamard gates are applied to selected qubits followed by 
#                       T gates. The state is converted to a density matrix and 
#                       then its Pauli expansion is computed.
#
#   :pbc        -> These are resource states arising from a reduction technique used in pbc.   
#                       A magic state is prepared by applying T gates (t_state).
#                       Then, for each hadamard pattern, Hadamard actions are applied
#                       directly on the T gate (i.e. on the magic state), and then 
#                       entanglement is generated via CZ operations.
#
# Results are saved to HDF5 files.
################################################################################

using Pkg
Pkg.activate(joinpath(dirname(@__DIR__), "MyProject"))
include(joinpath(@__DIR__, "robustness_functions.jl"))

using HDF5, Combinatorics

# ------------------------------------------------------------------------------
# Set up directory paths (relative to project root)
# ------------------------------------------------------------------------------
base_dir = dirname(@__DIR__)  # parent of local_robustness folder
data_folder = joinpath(base_dir, "data")
keys_folder = joinpath(base_dir, "keys")

# ------------------------------------------------------------------------------
# Parameters and Mode Selection
# ------------------------------------------------------------------------------

# Number of qubits
n = 3

# List of phase-space data file identifiers (adjust as needed)
file_names = ["deterministic_$n", "cnc_3_1"]

# List of phase-space data file identifiers (adjust as needed)
#file_names = ["deterministic_$n", "cnc_3_1"]

# Data folders (adjust paths as needed)
data_folder = "./data/"
output_folder = "./keys/"

# Mode selection:
#   :state_vector  -> Use state-vector approach.
#   :pauli         -> Use Pauli basis approach.
mode = :mbqc  # Change to :pauli to use the alternative approach

# ------------------------------------------------------------------------------
# Experiment Loop
# ------------------------------------------------------------------------------
for file_name in file_names
    # Construct HDF5 file path based on mode and qubit number.
    h5_file = joinpath(data_folder, file_name * (mode == :state_vector ? ".h5" : ".h5"))
    
    # Open and read the phase-space matrix.
    file = h5open(h5_file, "r")
    dataset_name = (mode == :state_vector ? file_name : file_name)
    R = read(file[dataset_name])
    close(file)
    
    println("Phase space: $file_name, n = $n")
    println("Loaded matrix size: ", size(R))
    
    if mode == :mbqc
        # State Vector Mode:
        # 1. Generate graph states from |+>^n.
        # 2. For each Hadamard pattern: apply H gates on selected qubits, then T gates.
        # 3. Convert the state to a density matrix and then to its Pauli expansion.
        graph_states, graph_names = generate_graph_state_vectors(n)
        for (idx, psi) in enumerate(graph_states)
            println("Graph: ", graph_names[idx])
            hadamard_bits = all_dit_strings(2, n)
            for bits in hadamard_bits
                state_transformed = copy(psi)
                hadamard_string = join(bits)
                println("Applying Hadamards: ", hadamard_string,"\n")
                for i in 1:n
                    if bits[i] == 1
                        state_transformed = hadamard_gate(n, i) * state_transformed
                    end
                    state_transformed = t_gate(n, i) * state_transformed
                end
                rho = state_transformed * (state_transformed')
                vector_rho = convert_matrix_to_pauli_basis(n, rho)
                exp_name = file_name * "_" * graph_names[idx] * "_h_" * hadamard_string * "_MBQC"
                println("Computing robustness for ", exp_name)
                robustness(R, vector_rho; threshold=1e-16, save=true, file_loc=keys_folder, file_name=exp_name)
                println("_________________________")
            end
            println("_________________________","\n")
        end
        
    elseif mode == :pbc
        # Pauli Mode:
        # For each Hadamard pattern:
        # 1. Start with a magic state (prepared via T gates).
        # 2. For each qubit with a 1 in the pattern, apply hadamard_action_on_pauli_basis.
        # 3. Then, generate graph states (via CZ operations) from the modified magic state.
        pauli_dict = PauliString(n)
        hadamard_bits = all_dit_strings(2, n)
        for bits in hadamard_bits
            magic_state = t_state(n)
            for i in 1:n
                if bits[i] == 1
                    magic_state = hadamard_action_on_pauli_basis(i, magic_state, pauli_dict)
                end
            end
            graph_states, graph_names = generate_graph_states(magic_state, n, pauli_dict)
            for (idx, state_transformed) in enumerate(graph_states)
                hadamard_string = join(bits)
                println("Graph: ", graph_names[idx], " with Hadamard pattern: ", hadamard_string,"\n")
                vector_rho = state_transformed
                exp_name = file_name * "_" * graph_names[idx] * "_h_" * hadamard_string * "_PBC"
                println("Computing robustness for ", exp_name)
                robustness(R, vector_rho; threshold=1e-16, save=true, file_loc=keys_folder, file_name=exp_name)
                println("_________________________")
            end
            println("_________________________")
        end
        
    else
        error("Invalid mode selected.")
    end
end
