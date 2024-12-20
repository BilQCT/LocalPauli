# initialize project:
using Pkg; Pkg.activate("./MyProject");;

# Load required packages
using HDF5, JuMP, GLPK, LinearAlgebra, Combinatorics

# Include utility functions
include("./libs/utils.jl"); include("./libs/pauli.jl")



function t_state(n)
    # Define single qubit t state
    t = Vector{Real}([1,1/sqrt(2),1/sqrt(2),0]);

    # initialize n-qubit t-state:
    T = Vector{Real}([1])
    for m in 1:n
        T = kron(T,t)
    end
    return Vector{Real}(T)
end

function cz_action_on_bsf(control,target,a::Vector{Int64})::Tuple{Vector{Int64},Int64}
    n = Int(length(a)/2)
    a_prime = copy(a)

    # extract control and target:
    x_control = a[control]; x_target = a[target]
    z_control = a[n+control]; z_target = a[n+target]

    # modify image:
    a_prime[n+control] = (z_control+x_target) % 2
    a_prime[n+target] = (z_target+x_control) % 2

    # phase flip:
    r =((x_control*x_target*z_control) + (x_control*x_target*z_target)) % 2
    return (a_prime,r)
end


function cz_action_on_pauli_basis(c,t,hermitian::Vector{Real},pauli_dict::PauliString)
    # Length of vector:
    N = length(hermitian)
    # check if dimension is 4^n: will raise error otherwise:
    n = Int(log2(N)/2)

    # initialize new hermitian:
    hermitian_prime = Vector{Real}([1 for i in 1:N])

    # for each i in 1:4^n:
    for i in 1:N
        # current pauli coefficient:
        pauli_coefficient = hermitian[i]
        # find current pauli:
        pauli = pauli_dict.int_to_bit[i]
        # find mapped pauli:
        pauli_prime, phase = cz_action_on_bsf(c,t,pauli)
        # mapped pauli coefficient:
        hermitian_prime[pauli_dict.bit_to_int[pauli_prime]] = pauli_coefficient*((-1)^phase)

        """
        # print results:
        if Bool(phase)
            println("$(pauli_dict.bit_to_pauli[pauli]) to -$(pauli_dict.bit_to_pauli[pauli_prime])")
        else
            println("$(pauli_dict.bit_to_pauli[pauli]) to $(pauli_dict.bit_to_pauli[pauli_prime])")
        end
        """
    end

    return hermitian_prime
end;;


function robustness(input_array, target_state; threshold=1e-16, save=true, file_loc="./", file_name="result")

    # Create a model with GLPK as the solver
    model = Model(GLPK.Optimizer)
    N = size(input_array)[2]

    # Define the decision variables
    @variable(model, x[1:N])
    @variable(model, u[1:N] >= 0)  # Auxiliary variables for the absolute values

    # Objective: Minimize the sum of u (which represents |x|)
    @objective(model, Min, sum(u))

    # Constraints for the absolute values
    @constraint(model, [i=1:N], u[i] >= x[i])
    @constraint(model, [i=1:N], u[i] >= -x[i])

    # Equality constraint: input_array * x = target_state
    @constraint(model, input_array * x .== target_state)

    # Solve the model
    optimize!(model)

    # Check the solver status
    if termination_status(model) != MOI.OPTIMAL
        error("Optimization failed with status: ", termination_status(model))
    end

    # Get the results
    x_value = value.(x)
    obj_value = objective_value(model)

    # Print the results
    println("Optimal objective value (min ||x||_1): ", obj_value)
    println("Optimal value of x: ", x_value)

    # Determine significant terms
    significant_terms = findall(x -> abs(x) > threshold, x_value)

    # Extract relevant terms
    decomposition_array = vcat(transpose(x_value[significant_terms]), input_array[:, significant_terms])

    if save
        # Save to HDF5 file
        h5_file = file_loc * file_name * ".h5"
        h5open(h5_file, "w") do file
            file[file_name] = decomposition_array
        end

        println("Optimal decomposition saved to $h5_file \n")
    end
end


#####################################################
#####################################################
# Compute robustness for 5-qubit magic cluster states:
#####################################################
#####################################################

using HDF5, Combinatorics

# Load the HDF5 file
file = h5open("./data/deterministic_5.h5", "r")

# Assign the dataset to R5
R5 = file["deterministic_5"][:,:]

# Close the file
close(file)

# Use the matrix
println("Loaded matrix size: ", size(R5))

# four qubit magic state:
n = 5
pauli_dict = PauliString(n)
T5 = t_state(n);


function generate_graph_states(magic_state, n, pauli_dict)
    # Initialize graph states
    Ln = copy(magic_state)
    for i in 1:n-1
        Ln = cz_action_on_pauli_basis(i, i+1, Ln, pauli_dict)
    end

    Cn = copy(magic_state)
    for i in 1:n
        if i < n
            Cn = cz_action_on_pauli_basis(i, i+1, Cn, pauli_dict)
        else
            Cn = cz_action_on_pauli_basis(1, n, Cn, pauli_dict)
        end
    end

    # Complete Graph
    pairs = [p for p in combinations(1:n, 2)]
    Kn = copy(magic_state)
    for p in pairs
        Kn = cz_action_on_pauli_basis(p[1], p[2], Kn, pauli_dict)
    end

    return [Ln, Cn, Kn],["L$n","C$n","K$n"]
end

# Call the function
graph_states, graph_names = generate_graph_states(T5, n, pauli_dict)

for i in 1:3
    # Identify graph:
    graph_state = graph_states[i]
    # Identify name:
    graph_name = graph_names[i]

    println("Computing robustness for $graph_name:\n")
    
    # Compute robustness:
    robustness(R5, graph_state; threshold=1e-16, save=true, file_loc="./", file_name="deterministic_5_"*graph_name)
end