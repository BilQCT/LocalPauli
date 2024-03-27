using LinearAlgebra
using Combinatorics
using YAML

include("stab.jl")

function do_commute(T_a::Vector{Int}, T_b::Vector{Int})
    """
    Check if two Pauli operators commute or not.

    Args:
    T_a: First Pauli operator as a binary symplectic vector.
    T_b: Second Pauli operator as a binary symplectic vector.

    Returns:
    bool: True if the two Pauli operators commute, False otherwise.
    """
    n = length(T_a) ÷ 2
    T_ax = T_a[1:n]
    T_az = T_a[n+1:end]

    T_bx = T_b[1:n]
    T_bz = T_b[n+1:end]

    return (dot(T_az, T_bx) + dot(T_ax, T_bz)) % 2 == 0
end

function get_pauli_string(T::Vector{Int})
    """
    Get the Pauli string representation of a Pauli operator.

    Args:
    T: Pauli operator as a binary symplectic vector.

    Returns:
    str: Pauli string representation of the Pauli operator.
    """
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

function find_isotropic_gens(isotropic::Set{Vector{Int}})::Set{Vector{Int}}
    """
    Find the isotropic generators of a given isotropic set of Pauli operators.

    Args:
    isotropic: Isotropic set of Pauli operators.

    Returns:
    Set{Vector{Int}}: Isotropic generators of the given isotropic set.
    """
     
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

mutable struct MaximalCncSet
    """
    Structure to store Closed maximal Non-contextual(CNC) sets of Pauli operators. By using the characterisation of 
    maximal CNC sets given in the Raussendorf et al, 'Phase space simulation method for quantum computation
    with magic states on qubits' paper, we can represent a maximal CNC set by storing the isotropic generators and
    anticommuting Paulis. 

    Attributes:
    isotropic_gens: Set of isotropic generators of the CNC set.
    anticommuting_paulis: Set of anticommuting Paulis of the CNC set. Here anticommuting paulis are a 
    representatives of the coset of I_a for each I_a in the CNC set.
    n: Number of qubits in the system.
    m: Parameter that determines the dimension of the isotropic set I. Dimension is n - m.
    """

    isotropic_gens::Set{Vector{Int}}
    anticommuting_paulis::Set{Vector{Int}}
    n::Int
    m::Int

    function MaximalCncSet(isotropic_gens::Set{Vector{Int}}, anticommuting_paulis::Set{Vector{Int}})
        """
        Constructor for the MaximalCncSet structure that creates a new maximal CNC set with the given isotropic
        generators and anticommuting Paulis.

        Args:
        isotropic_gens: Set of isotropic generators of the CNC set.
        anticommuting_paulis: Set of anticommuting Paulis of the CNC set.
        """
        n = length(first(isotropic_gens)) ÷ 2
        m = (length(anticommuting_paulis) - 1) ÷ 2
        new(isotropic_gens, anticommuting_paulis, n, m)
    end

    function MaximalCncSet(cnc_set::Set{Vector{Int}})
        """
        Constructor for the MaximalCncSet structure that creates a new maximal CNC set from a given the full set 
        of Pauli operators in the maximal CNC set. It extracts the isotropic generators and anticommuting Paulis
        of the full CNC set and stores them.

        Args:
        cnc_set: Full set of Pauli operators in the maximal CNC set.
        """
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

function Base.show(io::IO, cnc_set::MaximalCncSet)
    """
    Display the isotropic generators and anticommuting Paulis of a given maximal CNC set. In a basic sense, this
    changes what it will be displayed when the CNC set is printed.

    Args:
    io: Output stream.
    cnc_set: Maximal CNC set to display.
    """
    println(io, "Maximal CNC set:")
    println(io, "- Isotropic generators:")
    isotropic_gens_str = [get_pauli_string(p) for p in cnc_set.isotropic_gens]
    print(io, "  ")
    println(io, isotropic_gens_str)

    println(io, "- Anticommuting Paulis:")
    anticommuting_paulis_str = [get_pauli_string(p) for p in cnc_set.anticommuting_paulis]
    print(io, "  ")
    println(io, anticommuting_paulis_str)
end


mutable struct MaximalCnc
    """
    Structure to store a maximal CNC with a given Pauli operator set and value assignment. 

    Attributes:
    cnc_set: Maximal CNC set of the CNC.
    value_assignment: Value assignment of the Pauli operators in the CNC.
    """
    cnc_set::MaximalCncSet
    value_assignment::Dict{Vector{Int},Int}

    function MaximalCnc(cnc_set::MaximalCncSet, value_assignment::Dict{Vector{Int},Int})
        """
        Constructor for the MaximalCnc structure that creates a new maximal CNC with the given maximal CNC set
        and value assignment.
        """
        new(cnc_set, value_assignment)
    end
end

function generate_all_cncs_for_given_cnc_set(cnc_set::MaximalCncSet)::Set{MaximalCnc}
    """
    Generate all the possible maximal CNCs for a given maximal CNC set by finding all the possible value assignments.
    By characterisation the maximal CNCs and the properties of the value assignment, value assignment can be 
    determined by choosing the values of isotropic generators and anticommuting Paulis independently. For a maximal 
    CNC with (n,m), there are n - m generators for the isotropic subspace and 2m + 1 anticommuting Paulis. So, there
    are n-m + 2m + 1 = n + m + 1 values to assign to determine the value assignment. Hence there are 2^(n+m+1) 
    possible CNCs for a given maximal CNC set with (n,m).
    """

    # Create an array of all the independent Paulis in the maximal CNC set
    indep_paulis = collect(union(cnc_set.isotropic_gens, cnc_set.anticommuting_paulis))
    num_indep = length(indep_paulis)
    n = cnc_set.n

    cncs = Set{MaximalCnc}()
    
    # We need to assign each possible value for independent Paulis. For that, we can iterate over all binary strings
    # of length num_indep and determine the value assignment for independent Paulis using this array. For example,
    # if the binary string is 101, then the first and third independent Paulis will have a value of (-1)^1 = -1 
    # and the second independent Pauli will have a value of (-1)^0 = 1.
    # We are iterating from 0 to 2^num_indep-1 to get all the possible binary strings of length num_indep
    for i in 0:(2^num_indep-1)
        identity = [0 for i in 1:2n]

        # Initialize the value assignment by assigning the identity to 1
        value_assignment = Dict{Vector{Int},Int}(identity => 1)

        # Convert the integer to a binary string of length num_indep
        bitstring = string(i, base=2, pad=num_indep)
        value_array = [parse(Int, c) for c in bitstring]

        # Assign the values specified by binary array to the independent Paulis
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

        push!(cncs, MaximalCnc(cnc_set, value_assignment))
    end

    return cncs
end

function Base.show(io::IO, cnc::MaximalCnc)
    println(io, "CNC:")
    println(io, cnc.cnc_set)
    println(io, "- Value assignment:")
    for (p, v) in cnc.value_assignment
        println(io, "  $(get_pauli_string(p)): $v")
    end
end

function get_full_cnc_set(n::Int, m::Int)::Set{Vector{Int}}
    """
    Get the full set of Pauli operators for a canonical  maximal CNC set with (n,m) parameters.

    Args:
    n: Number of qubits in the system.
    m: Parameter that determines the dimension of the isotropic set I. Dimension is n - m.

    Returns:
    Set{Vector{Int}}: Full set of Pauli operators in the canonical maximal CNC set with (n,m) parameters.
    """
    data = YAML.load_file("cnc.yaml")
    return Set(data[n][m]["full_set"])
end