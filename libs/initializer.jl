using JuMP
const jmp = JuMP

using GLPK

include("cnc.jl")

function initial_lp(V::Vector{Float64},R)
    d = size(R)[1]; N = size(R)[2]

    # define model:
    model = jmp.Model(GLPK.Optimizer)

    # define variables: non-negative:
    @variable(model,x[1:N] >= 0)

    # define constraints:
    @constraint(model, cons, R * x .== V)

    # define objective function: We only want to find a solution:
    @objective(model, Min, 0)

    # solve model
    optimize!(model)


    if is_solved_and_feasible(model)
        mixture = "Convex"

        return model,x,mixture

    else
        mixture = "Affine"

        # define model:
        model = jmp.Model(GLPK.Optimizer)

        # define variables: non-negative:
        @variable(model,x[1:N] >= 0)

        # define constraints:
        @constraint(model, cons, R * x .== V)

        # define objective function: We only want to find a solution:
        @objective(model, Min, 0)

        # solve model
        optimize!(model)

        return model,x,mixture

    end
end


function approximate_input(V,e)
    W = [Float64(round(BigFloat(v);digits = e)) for v in V]
    return W
end


mutable struct Initializer
    """
    Structure to initialize simulation for Hermitian operator V given in terms of Pauli coefficients in lex order.

    InitVector: Initial vector in Pauli basis
    Feasible: Returns true if LP is feasible; i.e., V lies in convex hull of R
    Supports: CNC sets for which have positive coefficients in solution
    InitDist: Dictionary mapping CNC pairs (Omega,gamma) to their positive coefficients
    """
    InitVector::Vector{Float64}
    Feasible::Bool
    Mixture::String
    Supports::Vector{MaximalCnc}
    InitDist::Dict{MaximalCnc,Float64}

    function Initializer(n::Int,V::Vector{Float64},R)
        """
        Constructor for ...
        """
        PS = PauliString(n)

        InitVector = Vector{Float64}(V)

        # run model
        model, variables, mixture = initial_lp(V,R)
        Feasible = is_solved_and_feasible(model)
        Mixture = mixture
        
        # coefficients of convex mixture:
        values = value.(variables)
        nonzero = findall(x->x != 0,values)
        coefficients = values[nonzero]

        # vertices in mixture:
        Supports = [pauli_basis_to_cnc(R[:,i],PS) for i in nonzero]

        # Initial distribution:
        InitDist = Dict(zip(Supports,coefficients))

        new(InitVector,Feasible,Mixture,Supports,InitDist)
    end
end



