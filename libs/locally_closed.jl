using Combinatorics
using StatsBase

include("cnc.jl")

mutable struct AlmostMaximalTwoDimSet
    n::Int
    m::Int
    chosen_locals::Vector{Vector{Int}}
    independent_nonlocals::Vector{Vector{Int}}

    function AlmostMaximalTwoDimSet(n::Int, m::Int, chosen_locals::Vector{Vector{Int}}, independent_nonlocals::Vector{Vector{Int}})
        new(n, m, chosen_locals, independent_nonlocals)
    end
end

mutable struct AlmostMaximalTwoDim
    almost_maximal_two_dim_set::AlmostMaximalTwoDimSet
    value_assignment::Dict{Vector{Int},Int}

    function AlmostMaximalTwoDim(almost_maximal_two_dim_set::AlmostMaximalTwoDimSet, value_assignment::Dict{Vector{Int},Int})
        new(almost_maximal_two_dim_set, value_assignment)
    end
end

function generate_random_two_dim_almost_maximal_sets(n::Int, m::Int)::AlmostMaximalTwoDimSet
    """
    Generates a random almost maximal set of pauli operators in two dimensions with n local operators in Alice and 
    m local operators in Bob.

    Args:
        n: The number of local operators in Alice.
        m: The number of local operators in Bob.

    Returns:
        An almost maximal paulis set.
    """

    XI = [1, 0, 0, 0]
    YI = [1, 0, 1, 0]
    ZI = [0, 0, 1, 0]
    IX = [0, 1, 0, 0]
    IY = [0, 1, 0, 1]
    IZ = [0, 0, 0, 1]

    all_locals = [XI, YI, ZI, IX, IY, IZ]

    alice_locals = [XI, YI, ZI]
    bob_locals = [IX, IY, IZ]

    random_local_alice = sample(alice_locals, n, replace=false)
    random_local_bob = sample(bob_locals, m, replace=false)

    chosen_locals = vcat(random_local_alice, random_local_bob)

    other_locals = []
    for i in 1:6
        if !(all_locals[i] in chosen_locals)
            push!(other_locals, all_locals[i])
        end
    end

    independent_nonlocals::Vector{Vector{Int}} = []

    for used_locals in combinations(other_locals, 2)
        if do_commute(used_locals[1], used_locals[2])
            push!(independent_nonlocals, used_locals[1] + used_locals[2])
        end
    end

    return AlmostMaximalTwoDimSet(n, m, chosen_locals, independent_nonlocals)

end


function generate_all_two_dim_almost_maximals(almost_maximal::AlmostMaximalTwoDimSet)::Set{AlmostMaximalTwoDim}

    indep_paulis = vcat(almost_maximal.independent_nonlocals, almost_maximal.chosen_locals)
    num_indep = length(indep_paulis)

    almost_maximals = Set{AlmostMaximalTwoDim}()

    for i in 0:(2^num_indep-1)
        identity = [0, 0, 0, 0]

        # Initialize the value assignment by assigning the identity to 1
        value_assignment = Dict{Vector{Int},Int}(identity => 1)

        # Convert the integer to a binary string of length num_indep
        bitstring = string(i, base=2, pad=num_indep)
        value_array = [parse(Int, c) for c in bitstring]

        # Assign the values specified by binary array to the independent Paulis
        for (p, v) in zip(indep_paulis, value_array)
            value_assignment[p] = (-1)^v
        end

        for used_locals in combinations(almost_maximal.chosen_locals, 2)
            if do_commute(used_locals[1], used_locals[2])
                pauli = (used_locals[1] + used_locals[2]) .% 2
                value_assignment[pauli] = value_assignment[used_locals[1]] * value_assignment[used_locals[2]]
            end
        end

        push!(almost_maximals, AlmostMaximalTwoDim(almost_maximal, value_assignment))
    end

    return almost_maximals
end


function almost_maximal_to_pauli_basis(almost_maximal::AlmostMaximalTwoDim, ps::PauliStrings)
    bit_strings = ps.bit_strings
    N = length(bit_strings); V = []

    gamma = almost_maximal.value_assignment
    gamma_keys = collect(keys(gamma))

    for x in bit_strings
        if x in gamma_keys; push!(V,gamma[x]); else push!(V,0); end;
    end

    return V
end