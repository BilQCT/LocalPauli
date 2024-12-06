# initialize project:
using Pkg; Pkg.activate("./MyProject");;

# include symplectic orbit package:
include("./libs/symplectic.jl");; include("./libs/cnc.jl");;

n = 3;; m = 1;;

println("Generating maximal cnc sets for $(n) qubits with m=$(m) \n");;

println("Generating canonical maximal ($(n),$(m)) cnc set. \n")

# generators of canonical n,m cnc set:
canonical_cnc_n_m_gens = generate_canonical_cnc_with_n_m(n,m)

# full cnc set:
canonical_cnc_n_m = find_full_set_for_given_cnc_set(canonical_cnc_n_m_gens)

println("Generating symplectic group on $(n) qubits. \n")

# create symplectic group:
symplectic_n = SympPerm(n);;

println("Generating symplectic orbit of canonical maximal isotropic subspace. \n")

# Generate orbit
cnc_n_m_orbit = SympOrbit(n,symplectic_n,canonical_cnc_n_m).Orbit

println("Computation of orbit is completed. Now saving. \n")

# create output name
outputfile = "cnc_$(n)_$(m)_sets.jld"

using JLD2
@save outputfile cnc_n_m_orbit

println("File saved.")

