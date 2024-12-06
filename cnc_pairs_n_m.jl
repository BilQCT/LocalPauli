# initialize project:
using Pkg; Pkg.activate("./MyProject");;

# include symplectic orbit package:
include("./libs/symplectic.jl");; include("./libs/cnc.jl");;

n = 2;; m = 1;;

using JLD2

println("Loading maximal ($(n),$(m)) cnc sets.")

cnc_sets_n_m = collect(load("./data/cnc_$(n)_$(m)_sets.jld")["cnc_n_m_orbit"])

println("Determining generators of stabilizer and cosets. \n");;

cnc_gens_n_m = Vector{MaximalCncSet}([]);;
for cnc in cnc_sets_n_m
    push!(cnc_gens_n_m,MaximalCncSet(cnc))
end

println("Computation of generators is completed. Now saving. \n")

# create output name
outputfile = "cnc_gens_$(n)_$(m).yaml"

#using JLD2
#@save outputfile cnc_gens_n_m

using YAML

# Write the orbits to a YAML file
open(outputfile, "w") do file
    YAML.write(file, cnc_gens_n_m)
end

println("File saved.")

println("Determining value assignments for maximal ($(n),$(m)) cnc set. \n");;

cnc_pairs_n_m = Vector{Set{MaximalCnc}}([]);;
for cnc in cnc_gens_n_m
    push!(cnc_pairs_n_m,generate_all_cncs_for_given_cnc_set(cnc))
end

println("Computation of cnc pairs is completed. Now saving. \n")

# create output name
outputfile = "cnc_pairs_$(n)_$(m).yaml"

#using JLD2
#@save outputfile cnc_pairs_n_m

using YAML

# Write the orbits to a YAML file
open(outputfile, "w") do file
    YAML.write(file, cnc_pairs_n_m)
end

println("File saved.")

