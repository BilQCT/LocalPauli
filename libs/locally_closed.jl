using Combinatorics
using StatsBase
using LinearAlgebra

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

function generate_isotropic_from_gens(gens::Vector{Vector{Int}})::Set{Vector{Int}}
    """
    Generates the isotropic set from the generators.

    Args:
        gens: The generators of the isotropic set.

    Returns:
        The isotropic set.
    """

    n = length(gens[1]) / 2

    identity = [0 for i in 1:2n]
    isotropic = Set{Vector{Int}}([identity])

    for i in 1:length(gens)
        for used_gens in combinations(gens, i)
            pauli = identity
            for used_gen in used_gens
                pauli = (pauli + used_gen) .% 2
            end
            push!(isotropic, pauli)
        end
    end

    return isotropic
end

function generate_all_three_dim_maximal_local_isotropics_gens()::Set{Set{Vector{Int}}}

    XII = [1, 0, 0, 0, 0, 0]
    YII = [1, 0, 0, 1, 0, 0]
    ZII = [0, 0, 0, 1, 0, 0]

    alice_locals = [XII, YII, ZII]

    IXI = [0, 1, 0, 0, 0, 0]
    IYI = [0, 1, 0, 0, 1, 0]
    IZI = [0, 0, 0, 0, 1, 0]

    bob_locals = [IXI, IYI, IZI]

    IIX = [0, 0, 1, 0, 0, 0]
    IIY = [0, 0, 1, 0, 0, 1]
    IIZ = [0, 0, 0, 0, 0, 1]

    charlie_locals = [IIX, IIY, IIZ]

    all_maximal_local_isotropics_gens = Set{Set{Vector{Int}}}()

    for alice_local in alice_locals
        for bob_local in bob_locals
            for charlie_local in charlie_locals
                gens = Set{Vector{Int}}([alice_local, bob_local, charlie_local])
                push!(all_maximal_local_isotropics_gens, gens)
            end
        end
    end

    return all_maximal_local_isotropics_gens
end

function do_locally_commute(pauli1::Vector{Int}, pauli2::Vector{Int})::Bool
    
    if length(pauli1) != length(pauli2)
        return false
    end

    pauli_str_1 = get_pauli_string(pauli1)
    pauli_str_2 = get_pauli_string(pauli2)


    n = length(pauli_str_1)

    for i in 1:n
        if pauli_str_1[i] != pauli_str_2[i]
            if pauli_str_1[i] != 'I' && pauli_str_2[i] != 'I'
                return false
            end
        end
    end

    return true
end

function find_local_closure(omega::Set{Vector{Int}})
    is_closed = false

    while !is_closed
        is_closed = true
        for paulis in combinations(collect(omega), 2)
            pauli1 = paulis[1]
            pauli2 = paulis[2]
            if do_locally_commute(pauli1, pauli2) && !((pauli1 + pauli2) .%2 in omega)
                is_closed = false
                push!(omega, (pauli1 + pauli2) .% 2)
            end
        end
    end

    return omega
end

function is_locally_closed(omega::Set{Vector{Int}})
    if omega == find_local_closure(omega)
        return true
    end

    return false
end

function is_almost_maximal(omega::Set{Vector{Int}})
    l = length(first(omega))
    is_almost_maximal = true

    for i in 0:(2^l - 1)
        bitstring = string(i, base=2, pad=l)
        pauli = [parse(Int, c) for c in bitstring]

        if !(pauli in omega)
            new_omega = copy(omega)
            push!(new_omega, pauli)
            if length(find_local_closure(new_omega)) != 2^l
                println(get_pauli_string(pauli), " is added and closure is not maximal")
                return false
            end
        end
    end

    return true
end

function find_independent_paulis(omega::Set{Vector{Int}})
    if length(omega) == 0
        return []
    end
    omega = copy(omega)
    identity = [0 for i in 1:length(first(omega))]
    delete!(omega, identity)
    inferred_paulis = Set{Vector{Int}}([identity])

    random_indep_pauli = sample(collect(omega), 1)[1]
    delete!(omega, random_indep_pauli)
    push!(inferred_paulis, random_indep_pauli)
    independent_paulis = [random_indep_pauli]

    while length(omega) > 0
        random_indep_pauli = sample(collect(omega), 1)[1]
        delete!(omega, random_indep_pauli)
        newly_inferred_paulis = []
        for pauli in inferred_paulis
            if do_locally_commute(random_indep_pauli, pauli)
                inferred_pauli = (random_indep_pauli + pauli) .% 2
                push!(newly_inferred_paulis, inferred_pauli)
                delete!(omega, inferred_pauli)
            end
        end

        for pauli in newly_inferred_paulis
            push!(inferred_paulis, pauli)
        end

        push!(independent_paulis, random_indep_pauli)
    end

    return independent_paulis
end

function find_all_possible_local_value_assignments(omega::Set{Vector{Int}})
    n = length(first(omega)) / 2
    l = length(omega)

    omega = copy(omega)
    identity = [0 for i in 1:2n]
    delete!(omega, identity)

    all_value_assignments = Set{Dict{Vector{Int},Int}}()

    values_set = [i for i in 0:2^l-1]

    if l < 20
        values_set = [i for i in 0:2^l-1]
    else
        values_set = rand(0:2^l-1, 2^20)
    end

    for i in values_set
        # Initialize the value assignment by assigning the identity to 1
        value_assignment = Dict{Vector{Int},Int}(identity => 1)

        # Convert the integer to a binary string of length l
        bitstring = string(i, base=2, pad=l)
        value_array = [parse(Int, c) for c in bitstring]

        # Assign the values specified by binary array to the independent Paulis
        for (p, v) in zip(collect(omega), value_array)
            value_assignment[p] = (-1)^v
        end

        locally_value_assignment = true
        for paulis in combinations(collect(omega), 2)
            pauli1 = paulis[1]
            pauli2 = paulis[2]
            if do_locally_commute(pauli1, pauli2)
                pauli = (pauli1 + pauli2) .% 2
                if value_assignment[pauli] != value_assignment[pauli1] * value_assignment[pauli2]
                    locally_value_assignment = false
                    break
                end
            end
        end

        if locally_value_assignment
            push!(all_value_assignments, value_assignment)
        end
    end

    return all_value_assignments
end

function value_assignment_to_pauli_basis(value_assignment::Dict{Vector{Int},Int}, n::Int)
    ps = PauliStrings(n)
    bit_strings = ps.bit_strings
    V = []

    gamma_keys = collect(keys(value_assignment))

    for x in bit_strings
        if x in gamma_keys; push!(V,value_assignment[x]); else push!(V,0); end;
    end

    return V
end

function find_ranks_count_for_given_set_of_value_assignments(value_assignments::Set{Dict{Vector{Int},Int}}, n::Int, stab_coeffs)
    
    # initialize matrix:
    M = Array{Int}(undef, 4^n, 0)

    for value_assignment in value_assignments
        A = value_assignment_to_pauli_basis(value_assignment, n)
        M = hcat(M, A)
    end

    rank_counts = Dict{Int, Int}()
    H = stab_coeffs * M

    for j in 1:size(H)[2]
        Z = findall(x->x==0,H[:,j])
        AZ = stab_coeffs[Z,:]
        rank_AZ = rank(AZ)
        if haskey(rank_counts, rank_AZ)
            rank_counts[rank_AZ] += 1
        else
            rank_counts[rank_AZ] = 1
        end
    end

    return rank_counts
end