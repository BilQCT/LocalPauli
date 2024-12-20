# initialize project:
using Pkg; Pkg.activate("./MyProject");;
# Load required packages
using HDF5
# Include utility functions
include("./libs/utils.jl"); include("./libs/pauli.jl")

function generate_local_paulis(n::Int)
    pauli_operators = ["X", "Y", "Z"]
    n_qubit_local_paulis = String[]

    for pos in 1:n
        for op in pauli_operators
            # Create Pauli string with 'op' at 'pos' and 'I' elsewhere
            pauli_str = String(repeat('I', n))
            pauli_str = pauli_str[1:pos-1] * op * pauli_str[pos+1:end]
            push!(n_qubit_local_paulis, pauli_str)
        end
    end
    return n_qubit_local_paulis
end



function local_decomposition_pauli(n,pauli::Vector{Int})
    local_decomposition = Vector{String}([])
    for i in 1:n
        pauli_i = pauli[[i,i+n]]
        if pauli_i == Vector{Int}([1,0])
            A = "X"
            #push!(local_decomposition,"X")
        elseif pauli_i == Vector{Int}([0,1])
            A = "Z"
            #push!(local_decomposition,"Z")
        elseif pauli_i == Vector{Int}([1,1])
            A = "Y"
            #push!(local_decomposition,"Y")
        elseif pauli_i == Vector{Int}([0,0])
            A = "I"
            #push!(local_decomposition,"Y")
        end
        # Create Pauli string with 'op' at 'pos' and 'I' elsewhere
        pauli_str = String(repeat('I', n))
        pauli_str = pauli_str[1:i-1] * A * pauli_str[i+1:end]
        push!(local_decomposition, pauli_str)
    end
    return local_decomposition
end



function generate_all_deterministic_operators(n::Int64,pauli_dict::PauliString)
    # dimension
    dimension = 4^n
    # all possible outcomes:
    all_local_assignments = all_dit_strings(2,3*n)
    # local Paulis 
    local_paulis = generate_local_paulis(n)

    # initialize array:
    deterministic_array = Matrix{Int64}(undef,dimension,0)

    # generate vertex for every local assignment:
    for local_assignment in all_local_assignments
        # map local paulis to outcomes
        outcome_map = Dict(zip(local_paulis,local_assignment))
        outcome_map[String(repeat('I', n))] = 0
        # initialize vector:
        deterministic_vector = ones(Int64,dimension,1)

        # run through all paulis except identity:
        for i in 2:dimension
            pauli = pauli_dict.int_to_bit[i]
            local_decomposition_of_pauli = local_decomposition_pauli(n,pauli)
            values = [outcome_map[x] for x in local_decomposition_of_pauli]
            outcome = sum(values) % 2
            #println("Outcome for $(pauli_dict.int_to_pauli[i]): $(outcome)")
            if outcome == 1
                deterministic_vector[i] = -1
            end
        end
        deterministic_array = hcat(deterministic_array,deterministic_vector)
    end
    return deterministic_array
end




##########################################################
# Run script

# Check if a command-line argument is provided
if length(ARGS) == 0
    error("Please provide the number of qubits as a command-line argument.")
end

# Set the qubit number from the command line
n = parse(Int, ARGS[1])

# Create Pauli dictionary
pauli_dict = PauliString(n)

# Generate deterministic vertices
n_qubit_deterministic_operators = generate_all_deterministic_operators(n, pauli_dict)

# Save to HDF5 file
h5_file = "./data/deterministic_$(n).h5"
h5open(h5_file, "w") do file
    file["deterministic_$(n)"] = n_qubit_deterministic_operators
end

println("Matrix saved to $h5_file")


