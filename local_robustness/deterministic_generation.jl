using Pkg

# Activate the project environment located in "MyProject" folder
base_dir = dirname(@__DIR__)            # one level up from local_robustness
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using HDF5

# Set up directory paths relative to the project root
data_folder = joinpath(base_dir, "data")


# Generate deterministic vertices:
global D1 = Matrix{Int8}(undef,4,0)

# Find all possible 3-bit bitstrings:
combinations = collect(Iterators.product(0:1, 0:1, 0:1))

# Generate all vectors:
for combo in combinations
    local D = Int8[1, (-1)^combo[1], (-1)^combo[2], (-1)^combo[3]]
    global D1 = hcat(D1,D)
end


function generate_all_deterministic_operators(n::Int8)
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


# Explicit bounded default; the previous hard-coded n=7 allocation was enormous.
if abspath(PROGRAM_FILE) == @__FILE__
    for n in (isempty(ARGS) ? (2:4) : parse.(Int, ARGS))
        1 <= n <= 5 || error("Dense deterministic generation is restricted to n <= 5")
        data = generate_all_deterministic_operators(Int8(n))
        h5open(joinpath(data_folder, "deterministic_$(n).h5"), "w") do file
            file["data"] = data
            attributes(file)["pauli_convention"] = "lexicographic I,X,Y,Z; Tr(P A)"
        end
    end
end
