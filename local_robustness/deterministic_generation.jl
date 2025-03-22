using Pkg

# Activate the project environment located in "MyProject" folder
base_dir = dirname(@__DIR__)            # one level up from local_robustness
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5

# Set up directory paths relative to the project root
data_folder = joinpath(base_dir, "data")


# Generate deterministic vertices:
global D1 = Matrix{Int64}(undef,4,0)

# Find all possible 3-bit bitstrings:
combinations = collect(Iterators.product(0:1, 0:1, 0:1))

# Generate all vectors:
for combo in combinations
    local D = vcat([1],[(-1)^a for a in combo])
    global D1 = hcat(D1,D)
end


function generate_all_deterministic_operators(n::Int64)
    global D = copy(D1)
    if n == 1
        return D
    else
        for _ in 2:n
            D = kron(D,D1)
        end
    end
    return D
end


# qubit number
n = 1

println("Generating deterministic vertices.\n")

n_qubit_deterministic_operators = generate_all_deterministic_operators(n)

println("Saving deterministic vertices.\n")

# Save to HDF5 file
file_name = "deterministic_$(n)"
h5_file = joinpath(data_folder, file_name*".h5")
h5open(h5_file, "w") do file
    file[file_name] = n_qubit_deterministic_operators
end

println("Matrix saved to $h5_file\n")




