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