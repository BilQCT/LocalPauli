using LinearAlgebra

using Combinatorics

include("stab.jl")

function do_commute(T_a::Vector{Int}, T_b::Vector{Int})
    n = length(T_a) ÷ 2
    T_ax = T_a[1:n]
    T_az = T_a[n+1:end]

    T_bx = T_b[1:n]
    T_bz = T_b[n+1:end]

    return (dot(T_az, T_bx) + dot(T_ax, T_bz)) % 2 == 0
end

function get_pauli_string(T::Vector{Int})
    n = length(T) ÷ 2
    pauli_str = ""

    for i in 1:n
        if T[i] == 1 && T[i+n] == 1
            pauli_str *= "Y"
        elseif T[i+n] == 1
            pauli_str *= "Z"
        elseif T[i] == 1
            pauli_str *= "X"
        else
            pauli_str *= "I"
        end
    end

    return pauli_str
end

function find_isotropic_gens(isotropic::Set{Vector{Int}})
    n = length(first(isotropic)) ÷ 2
    isotropic = copy(isotropic)
    num_gens = round(Int, log2(length(isotropic)))
    generated_set = Set{Vector{Int}}([[0 for i in 1:2n]])
    newly_generated_set = Set{Vector{Int}}([[0 for i in 1:2n]])
    isotropic_gens = Set{Vector{Int}}()

    count = 0
    while count < num_gens
        for p in newly_generated_set
            delete!(isotropic, p)
        end
        if length(isotropic) == 0
            break
        end

        b = pop!(isotropic)
        push!(isotropic_gens, b)
        count += 1
        newly_gen_list = [(a + b) .% 2 for a in generated_set]
        gen_list = append!(newly_gen_list, [a for a in generated_set])
        newly_generated_set = Set{Vector{Int}}(newly_gen_list)
        generated_set = Set{Vector{Int}}(gen_list)
    end

    return isotropic_gens
end

mutable struct CncSet
    isotropic_gens::Set{Vector{Int}}
    anticommuting_paulis::Set{Vector{Int}}
    n::Int
    m::Int

    function CncSet(isotropic_gens::Set{Vector{Int}}, anticommuting_paulis::Set{Vector{Int}})
        n = length(first(isotropic_gens)) ÷ 2
        m = (length(anticommuting_paulis) - 1) ÷ 2
        new(isotropic_gens, anticommuting_paulis, n, m)
    end

    function CncSet(cnc_set::Set{Vector{Int}})
        isotropic = Set{Vector{Int}}()
        n = length(first(cnc_set)) ÷ 2

        for p in cnc_set
            if all(p2 -> do_commute(p, p2), cnc_set)
                push!(isotropic, p)
            end
        end

        for p in isotropic
            delete!(cnc_set, p)
        end

        if length(cnc_set) == 0
            anticommuting_paulis = Set{Vector{Int}}()
        else
            a = pop!(cnc_set)
            anticommuting_paulis = Set{Vector{Int}}([a])
            while length(cnc_set) > 0
                for p in anticommuting_paulis
                    for q in cnc_set
                        if do_commute(p, q)
                            delete!(cnc_set, q)
                        end
                    end
                end
                if length(cnc_set) == 0
                    break
                end
                b = pop!(cnc_set)
                push!(anticommuting_paulis, b)
            end
        end

        isotropic_gens = find_isotropic_gens(isotropic)

        if length(isotropic_gens) == 0
            isotropic_gens = Set{Vector{Int}}()
        end

        m = (length(anticommuting_paulis) - 1) ÷ 2
        new(isotropic_gens, anticommuting_paulis, n, m)
    end
end

function Base.show(io::IO, cnc_set::CncSet)
    println(io, "CNC set:")
    println(io, "- Isotropic generators:")
    isotropic_gens_str = [get_pauli_string(p) for p in cnc_set.isotropic_gens]
    print("  ")
    println(io, isotropic_gens_str)

    println(io, "- Anticommuting Paulis:")
    anticommuting_paulis_str = [get_pauli_string(p) for p in cnc_set.anticommuting_paulis]
    print("  ")
    println(io, anticommuting_paulis_str)
end


mutable struct CNC
    cnc_set::CncSet
    value_assignment::Dict{Vector{Int},Int}

    function CNC(cnc_set::CncSet, value_assignment::Dict{Vector{Int},Int})
        new(cnc_set, value_assignment)
    end
end

function generate_all_cncs_for_given_cnc_set(cnc_set::CncSet)::Set{CNC}
    # We can choose the value of isotropic generators and anticommuting Paulis independently in each value assignment
    # Note that anticommuting_paulis stores only the a for each I_a in the CNC set
    num_indep = length(cnc_set.isotropic_gens) + length(cnc_set.anticommuting_paulis)

    indep_paulis = collect(union(cnc_set.isotropic_gens, cnc_set.anticommuting_paulis))
    n = cnc_set.n

    cncs = Set{CNC}()
    for i in 0:(2^num_indep-1)
        identity = [0 for i in 1:2n]
        value_assignment = Dict{Vector{Int},Int}(identity => 1)

        # Generate the values of mentioned paulis in the value assignment as a value_array
        bitstring = string(i, base=2, pad=num_indep)
        value_array = [parse(Int, c) for c in bitstring]

        # Assign the values to the isotropic generators and anticommuting Paulis
        for (p, v) in zip(indep_paulis, value_array)
            value_assignment[p] = (-1)^v
        end

        non_identity_isotropic = Set{Vector{Int}}(cnc_set.isotropic_gens)
        # Extend the value assignment to all the Paulis in the isotropic
        for i in 2:(length(cnc_set.isotropic_gens))
            for used_generators in combinations(collect(cnc_set.isotropic_gens), i)
                value = 1
                pauli = [0 for i in 1:2n]
                for gen in used_generators
                    value *= value_assignment[gen]
                    pauli = (pauli + gen) .% 2
                end
                value_assignment[pauli] = value
                push!(non_identity_isotropic, pauli)
            end
        end

        # Extend the value assignment to all the Paulis in the anticommuting sets (I_a)
        for anti_com in cnc_set.anticommuting_paulis
            for p in non_identity_isotropic
                pauli = (p + anti_com) .% 2
                value_assignment[pauli] = value_assignment[p] * value_assignment[anti_com]
            end
        end

        push!(cncs, CNC(cnc_set, value_assignment))
    end

    return cncs
end

function Base.show(io::IO, cnc::CNC)
    println(io, "CNC:")
    println(io, cnc.cnc_set)
    println(io, "- Value assignment:")
    for (p, v) in cnc.value_assignment
        println(io, "  $(get_pauli_string(p)): $v")
    end
end

cnc = CncSet(Set{Vector{Int}}([[0, 1, 0, 0]]), Set{Vector{Int}}([[1, 0, 0, 0], [1, 0, 1, 0], [0, 0, 1, 0]]))

cncs = generate_all_cncs_for_given_cnc_set(cnc)
print("\n\n")
println(cncs)