# initialize project:
using Pkg; Pkg.activate("./MyProject");;

# include symplectic orbit package:
include("./libs/symplectic.jl"); include("./libs/cnc.jl");;

n = 4;; m_values = [1];;


println("Generating maximal cnc sets for $(n) qubits with m=$(m_values[1]) \n");;

cnc_n_m_values = generate_all_cncs(n, m_values);;

println("Computation of orbit is completed. Now saving. \n")



# create output name
outputfile = "cnc_$(n)_$(m_values[1])";;

using JLD2
@save outputfile cnc_n_m_values

println("File saved.")
