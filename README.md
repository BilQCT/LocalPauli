# StabTheory
Julia code for generating Stabilizer states and Clifford group

clifford_group.g: GAP script for generating clifford group and/or symplectic group. The bilinear form for Sp(2n,2) is modified to correspond to the quantum computing convention, rather than default GAP bilinear form. There is also code for induced group action of the Clifford group on the Pauli coefficients of a Hermitian matrix. 

clifford_group.jl: julia script for calling GAP and running GAP clifford_group.g functions.

stab.jl: julia script for generating stabilizer states.



# Requirements:

GAP
Polymake
Nemo
Combinatorics
YAML
JuMP
GLPK
Distributions

## Instructions

To create an environment for running StabTheory, in Julia REPL use the following commands:

using Pkg
Pkg.generate("path/to/MyProject")
Pkg.activate("path/to/MyProject")

Pkg.add("YourPackage")
