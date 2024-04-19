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

III = [0, 0, 0, 0, 0, 0]
XII = [1, 0, 0, 0, 0, 0]
YII = [1, 0, 0, 1, 0, 0]
ZII = [0, 0, 0, 1, 0, 0]
IXI = [0, 1, 0, 0, 0, 0]
IYI = [0, 1, 0, 0, 1, 0]
IZI = [0, 0, 0, 0, 1, 0]
IIX = [0, 0, 1, 0, 0, 0]
IIY = [0, 0, 1, 0, 0, 1]
IIZ = [0, 0, 0, 0, 0, 1]

IXX = (IXI + IIX) .% 2
IXY = (IXI + IIY) .% 2
IXZ = (IXI + IIZ) .% 2
IYX = (IYI + IIX) .% 2
IYY = (IYI + IIY) .% 2
IYZ = (IYI + IIZ) .% 2
IZX = (IZI + IIX) .% 2
IZY = (IZI + IIY) .% 2
IZZ = (IZI + IIZ) .% 2

XIX = (XII + IIX) .% 2
XIY = (XII + IIY) .% 2
XIZ = (XII + IIZ) .% 2
XXI = (XII + IXI) .% 2
XXX = (XXI + IIX) .% 2
XXY = (XXI + IIY) .% 2
XXZ = (XXI + IIZ) .% 2
XYI = (XII + IYI) .% 2
XYX = (XYI + IIX) .% 2
XYY = (XYI + IIY) .% 2
XYZ = (XYI + IIZ) .% 2
XZI = (XII + IZI) .% 2
XZX = (XZI + IIX) .% 2
XZY = (XZI + IIY) .% 2
XZZ = (XZI + IIZ) .% 2

YIX = (YII + IIX) .% 2
YIY = (YII + IIY) .% 2
YIZ = (YII + IIZ) .% 2
YXI = (YII + IXI) .% 2
YXX = (YXI + IIX) .% 2
YXY = (YXI + IIY) .% 2
YXZ = (YXI + IIZ) .% 2
YYI = (YII + IYI) .% 2
YYX = (YYI + IIX) .% 2
YYY = (YYI + IIY) .% 2
YYZ = (YYI + IIZ) .% 2
YZI = (YII + IZI) .% 2
YZX = (YZI + IIX) .% 2
YZY = (YZI + IIY) .% 2
YZZ = (YZI + IIZ) .% 2

ZIX = (ZII + IIX) .% 2
ZIY = (ZII + IIY) .% 2
ZIZ = (ZII + IIZ) .% 2
ZXI = (ZII + IXI) .% 2
ZXX = (ZXI + IIX) .% 2
ZXY = (ZXI + IIY) .% 2
ZXZ = (ZXI + IIZ) .% 2
ZYI = (ZII + IYI) .% 2
ZYX = (ZYI + IIX) .% 2
ZYY = (ZYI + IIY) .% 2
ZYZ = (ZYI + IIZ) .% 2
ZZI = (ZII + IZI) .% 2
ZZX = (ZZI + IIX) .% 2
ZZY = (ZZI + IIY) .% 2
ZZZ = (ZZI + IIZ) .% 2