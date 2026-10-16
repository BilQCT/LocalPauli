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
    bits = Int.(collect(x))
    iseven(length(bits)) || throw(ArgumentError("Pauli bit string must have even length"))
    all(b -> b in (0, 1), bits) || throw(ArgumentError("Pauli entries must be binary"))
    n = length(bits) ÷ 2
    local_paulis = (ComplexF64[1 0; 0 1], ComplexF64[0 1; 1 0],
                    ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
    result = ones(ComplexF64, 1, 1)
    for q in 1:n
        a, b = bits[q], bits[n+q]
        label = a == 0 ? (b == 0 ? 1 : 4) : (b == 0 ? 2 : 3)
        result = kron(result, local_paulis[label])
    end
    return result
end



############################################################################################
############################################################################################
############################################################################################



mutable struct PauliString
    bit_strings::Vector{Vector{Int}}
    pauli_strings::Vector{String}
    pauli_to_bit::Dict{String, Vector{Int}}
    bit_to_pauli::Dict{Vector{Int}, String}
    bit_to_int::Dict{Vector{Int}, Int}
    int_to_bit::Dict{Int, Vector{Int}}
    pauli_to_int::Dict{String, Int}
    int_to_pauli::Dict{Int, String}

    function PauliString(n::Int)
        bit_strings = vec(all_dit_strings(2, 2*n))
        bit_strings = [Int.(collect(x)) for x in bit_strings]

        pauli_strings = [get_pauli_string(x) for x in bit_strings]

        perm = sortperm(pauli_strings)
        pauli_strings = pauli_strings[perm]
        bit_strings = bit_strings[perm]

        pauli_to_bit = Dict(zip(pauli_strings, bit_strings))
        bit_to_pauli = Dict(zip(bit_strings, pauli_strings))

        bit_to_int = Dict(zip(bit_strings, 1:length(bit_strings)))
        int_to_bit = Dict(zip(1:length(bit_strings), bit_strings))

        pauli_to_int = Dict(
            p => bit_to_int[pauli_to_bit[p]]
            for p in pauli_strings
        )

        int_to_pauli = Dict(
            i => bit_to_pauli[int_to_bit[i]]
            for i in 1:length(bit_strings)
        )

        new(
            bit_strings,
            pauli_strings,
            pauli_to_bit,
            bit_to_pauli,
            bit_to_int,
            int_to_bit,
            pauli_to_int,
            int_to_pauli,
        )
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
