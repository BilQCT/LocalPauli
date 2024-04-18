# PauliString struct and related functions:

"""
Generate all dit strings recursively.

Parameters:
- `d`: The number of outcomes.
- `N`: The number of generators.
- `n`: The current iteration.
- `s`: The current set of dit strings.

Returns:
- All dit strings.
"""
function generate_dit_strings(d, N, n, s)
    if n < N
        Zd = [(i - 1) for i in 1:d]
        dit_strings = []
        for i in 1:length(s)
            for x in Zd
                push!(dit_strings, push!(copy(s[i]), x))
            end
        end
        s = dit_strings
        n = n + 1
        generate_dit_strings(d, N, n, s)
    else
        return s
    end
end


"""
Generate all possible dit strings.

Parameters:
- `d`: The number of outcomes.
- `N`: The number of generators.

Returns:
- All dit strings.
"""
function all_dit_strings(d, N)
    s = [[i - 1] for i in 1:d]
    n = 1
    return generate_dit_strings(d, N, n, s)
end

mutable struct PauliString
    """
    Structure to store a pauli strings, bit strings, and dictionary mapping between the two 

    Attributes:
    pauli_strings: for each operator (e.g., XI) a string 'xi'
    bit_strings: for each operator (e.g., XI) a bit string [1,0,0,0]
    pauli_to_bit: dictionary from pauli strings to bit strings
    bit_to_pauli: dictionary from bit strings to pauli strings
    bit_to_int: dictionary from bit string to integer order of pauli strings in lexicographic order
    """
    bit_strings::Vector{Vector{Int}}
    pauli_strings::Vector{String}
    pauli_to_bit::Dict{String,Vector{Int}}
    bit_to_pauli::Dict{Vector{Int},String}
    bit_to_int::Dict{Vector{Int},Int}
    int_to_bit::Dict{Int,Vector{Int}}

    function PauliString(n::Int)
        """
        Constructor for the PauliString structure that creates a vectors of pauli strings and bit strings and dictionaries mapping between them
        """
        bit_strings = all_dit_strings(2,2*n)
        pauli_strings = [get_pauli_string(x) for x in bit_strings]
        
        # In lex order
        perm = sortperm(pauli_strings)
        pauli_strings = pauli_strings[perm]; bit_strings = bit_strings[perm];

        pauli_to_bit = Dict(zip(pauli_strings,bit_strings))
        bit_to_pauli = Dict(zip(bit_strings,pauli_strings))
        bit_to_int = Dict(zip(bit_strings,[i for i in 1:length(bit_strings)]))
        int_to_bit = Dict(zip([i for i in 1:length(bit_strings)],bit_strings))
        new(bit_strings,pauli_strings,pauli_to_bit,bit_to_pauli,bit_to_int,int_to_bit)
    end
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