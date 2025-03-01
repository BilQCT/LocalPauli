############################################################################################
############################################################################################
############################################################################################

"""
Functions related to Pauli operators, Pauli group. 
"""


############################################################################################
############################################################################################
############################################################################################

function tensor_prod(A, n)
    """
    Compute the tensor product of a matrix `A` with itself `n` times.

    Parameters:
    - `A`: The matrix to be tensor producted.
    - `n`: The number of times to tensor product `A` with itself.

    Returns:
    - The resulting tensor product matrix.
    """
    An = A
    for i in 1:(n - 1)
        An = kron(An, A)
    end
    return An
end

function pauli(x)
    """
    Compute the Pauli operator corresponding to a given bit string.

    Parameters:
    - `x`: The bit string.

    Returns:
    - The Pauli operator corresponding to the input bit string.
    """
    if typeof(length(x) // 2) != Rational{Int64}
        return "Must be a bit-string of even length."
    elseif typeof(length(x) // 2) == Rational{Int64}
        n = Int(length(x) // 2)
        X = Matrix([0 1; 1 0])
        Z = Matrix([1 0; 0 -1])
        a = x[1:n]
        b = x[n + 1:end]
        if n == 1
            XaZb = (X^a[1] * Z^b[1])
        else
            XaZb = (X^a[1] * Z^b[1])
            for i in 2:n
                XaZb = kron(XaZb, (X^a[i] * Z^b[i]))
            end
        end
        arg_phase = (transpose(a) * b) % 4
        phase = im ^ (arg_phase)
        return phase * XaZb
    end
end



############################################################################################
############################################################################################
############################################################################################



# PauliString struct and related functions:

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
    pauli_to_int::Dict{String,Int}
    int_to_pauli::Dict{Int,String}

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
        pauli_to_int = Dict(zip(pauli_strings,[bit_to_int[pauli_to_bit[x]] for x in pauli_strings]))
        into_to_pauli = Dict(zip([i for i in 1:length(bit_strings)],[bit_to_pauli[int_to_bit[i]] for i in 1:length(bit_strings)]))
        new(bit_strings,pauli_strings,pauli_to_bit,bit_to_pauli,bit_to_int,int_to_bit,pauli_to_int,into_to_pauli)
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
