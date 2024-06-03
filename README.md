# StabTheory
Julia code for generating Stabilizer states and Clifford group

clifford_group.g: GAP script for generating clifford group and/or symplectic group. The bilinear form for Sp(2n,2) is modified to correspond to the quantum computing convention, rather than default GAP bilinear form. There is also code for induced group action of the Clifford group on the Pauli coefficients of a Hermitian matrix. 

clifford_group.jl: julia script for calling GAP and running GAP clifford_group.g functions.

stab.jl: julia script for generating stabilizer states.



# Requirements:

name = "MyProject"
uuid = "0758bb20-c2cc-4080-aeb8-35ec505a5029"
authors = ["selmanipek <selman.ipek.us@gmail.com>"]
version = "0.1.0"

[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
GAP = "c863536a-3901-11e9-33e7-d5cd0df7b904"
GraphRecipes = "bd48cda9-67a9-57be-86fa-5b3c104eda73"
Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Polymake = "d720cf60-89b5-51f5-aff5-213f193123e7"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
YAML = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"

## Instructions

To create an environment for running LocalLambda, in Julia REPL use the following commands:

using Pkg
Pkg.generate("path/to/MyProject")
Pkg.activate("path/to/MyProject")

Pkg.add("YourPackage")
