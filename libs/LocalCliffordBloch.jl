module LocalCliffordBloch

using GAP
using LinearAlgebra
using Combinatorics
include("MaxWeightSeed.jl")
using .MaxWeightSeed: gamma_elementary_p, seed_coordinates

export PauliString,
       pauli_n,
       Pauli_group,
       LocalCliffordGroup,
       local_clifford_orbit_of_point,
       local_clifford_gate_action_on_point,
       max_weight_bloch_vector,
       all_max_weight_fiducials,
       z3_to_local_symplectic,
       gamma_elementary_p

const g = GAP.Globals
const gjl = GAP.gap_to_julia
const gobj = GAP.GapObj

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

function all_dit_strings(d::Int, n::Int)
    return [collect(t) for t in Iterators.product(ntuple(_ -> 0:d-1, n)...)]
end

function z3_to_local_symplectic(u::Vector{Int})
    all(x -> 0 <= x <= 2, u) || throw(ArgumentError("inputs must lie in 0:2"))
    return vcat(Int.(u .!= 2), Int.(u .!= 0))
end

function max_weight_bloch_vector(n::Int, p::Int; pauli_info=PauliString(n))
    return seed_coordinates(n, p)
end

function all_max_weight_fiducials(n::Int; pauli_info=PauliString(n))
    return Dict(p => max_weight_bloch_vector(n, p; pauli_info=pauli_info) for p in 2:n)
end

reverse_map(D::Dict) = Dict(v => k for (k, v) in D)

function _as_int_vector(x)
    if x isa GAP.GapObj
        return Int.(collect(gjl(x)))
    else
        return Int.(collect(x))
    end
end

# -----------------------------------------------------------------------------
# Pauli strings and matrices
# Convention: binary symplectic vector [a_1,...,a_n,b_1,...,b_n]
# represents T_(a,b) = i^(a⋅b) X^a Z^b.
# -----------------------------------------------------------------------------

function get_pauli_string(T::Vector{Int})
    n = length(T) ÷ 2
    pauli_str = ""

    for i in 1:n
        if T[i] == 1 && T[n+i] == 1
            pauli_str *= "Y"
        elseif T[n+i] == 1
            pauli_str *= "Z"
        elseif T[i] == 1
            pauli_str *= "X"
        else
            pauli_str *= "I"
        end
    end

    return pauli_str
end

mutable struct PauliString
    bit_strings::Vector{Vector{Int}}
    pauli_strings::Vector{String}
    pauli_to_bit::Dict{String,Vector{Int}}
    bit_to_pauli::Dict{Vector{Int},String}
    bit_to_int::Dict{Vector{Int},Int}
    int_to_bit::Dict{Int,Vector{Int}}
    pauli_to_int::Dict{String,Int}
    int_to_pauli::Dict{Int,String}
end

function PauliString(n::Int)
    # All binary symplectic vectors in E_n = Z_2^{2n}
    bit_strings = vec(all_dit_strings(2, 2*n))
    bit_strings = [Int.(collect(a)) for a in bit_strings]

    pauli_strings = [get_pauli_string(a) for a in bit_strings]

    # Lexicographic Pauli order: I < X < Y < Z, qubit 1 first.
    perm = sortperm(pauli_strings)

    pauli_strings = pauli_strings[perm]
    bit_strings = bit_strings[perm]

    pauli_to_bit = Dict(zip(pauli_strings, bit_strings))
    bit_to_pauli = Dict(zip(bit_strings, pauli_strings))
    bit_to_int = Dict(zip(bit_strings, 1:length(bit_strings)))
    int_to_bit = Dict(zip(1:length(bit_strings), bit_strings))
    pauli_to_int = Dict(p => bit_to_int[pauli_to_bit[p]] for p in pauli_strings)
    int_to_pauli = Dict(i => bit_to_pauli[int_to_bit[i]] for i in 1:length(bit_strings))

    return PauliString(
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

function pauli_n(ab)::GAP.GapObj
    abv = _as_int_vector(ab)

    if !iseven(length(abv))
        error("Pauli bit string must have even length.")
    end

    n = length(abv) ÷ 2
    a = abv[1:n]
    b = abv[n+1:2*n]

    X = gobj([0 1; 1 0])
    Z = g.DiagonalMat(gobj([1, -1]))

    M = X^a[1] * Z^b[1]
    if n > 1
        for i in 2:n
            M = g.KroneckerProduct(M, X^a[i] * Z^b[i])
        end
    end

    phase = g.E(4)^mod(sum(a .* b), 4)
    return phase * M
end

function Pauli_group(n::Int)
    gens = GapObj[]

    for i in 1:n
        x = zeros(Int, 2*n)
        x[i] = 1

        z = zeros(Int, 2*n)
        z[n+i] = 1

        y = mod.(x + z, 2)

        push!(gens, pauli_n(gobj(x)))
        push!(gens, pauli_n(gobj(y)))
        push!(gens, pauli_n(gobj(z)))
    end

    return g.Group(gobj(gens))
end

# -----------------------------------------------------------------------------
# Basic Clifford gates over GAP cyclotomics
# -----------------------------------------------------------------------------

const R2 = g.E(8) + g.E(8)^7
const Id2 = gobj([1 0; 0 1])
const Phase = gobj([1 0; 0 g.E(4)])
const Had = (1/R2) * gobj([1 1; 1 -1])

function n_fold_tensor(n::Int, mat)
    if n == 0
        return g.DiagonalMat(gobj([1]))
    elseif n == 1
        return mat
    elseif n > 1
        result = mat
        for _ in 1:n-1
            result = g.KroneckerProduct(result, mat)
        end
        return result
    else
        error("n must be nonnegative.")
    end
end

function hadamard_j(j::Int, n::Int)
    Hj = (j == 1) ? Had : Id2
    for i in 2:n
        Hj = g.KroneckerProduct(Hj, i == j ? Had : Id2)
    end
    return Hj
end

function phase_j(j::Int, n::Int)
    Pj = (j == 1) ? Phase : Id2
    for i in 2:n
        Pj = g.KroneckerProduct(Pj, i == j ? Phase : Id2)
    end
    return Pj
end

function cnot(c::Int, t::Int, n::Int)
    if n == 1
        return Id2
    end

    Proj00 = g.DiagonalMat(gobj([1, 0]))
    Proj11 = g.DiagonalMat(gobj([0, 1]))
    NOT = gobj([0 1; 1 0])

    op0 = GAP.GapObj[]
    op1 = GAP.GapObj[]

    for i in 1:n
        if i == c
            push!(op0, Proj00)
            push!(op1, Proj11)
        elseif i == t
            push!(op0, Id2)
            push!(op1, NOT)
        else
            push!(op0, Id2)
            push!(op1, Id2)
        end
    end

    A = op0[1]
    B = op1[1]
    for i in 2:n
        A = g.KroneckerProduct(A, op0[i])
        B = g.KroneckerProduct(B, op1[i])
    end

    return A + B
end

function swap(c::Int, t::Int, n::Int)
    hc_ht = hadamard_j(c, n) * hadamard_j(t, n)
    if c < t
        cnot_ct = cnot(c, t, n)
        cnot_tc = hc_ht * cnot_ct * hc_ht
        return cnot_ct * cnot_tc * cnot_ct
    else
        cnot_tc = cnot(t, c, n)
        cnot_ct = hc_ht * cnot_tc * hc_ht
        return cnot_tc * cnot_ct * cnot_tc
    end
end

function local_clifford_unitaries(n::Int)::GAP.GapObj
    gens = GAP.GapObj[]

    for i in 1:n
        push!(gens, hadamard_j(i, n))
        push!(gens, phase_j(i, n))
        if i < n
            push!(gens, swap(i, i + 1, n))
        end
    end

    return g.Group(gobj(gens))
end

# -----------------------------------------------------------------------------
# GAP group wrappers and local Clifford action on signed Paulis
# -----------------------------------------------------------------------------

mutable struct GapGroupInfo
    Grp::GAP.GapObj
    Order::Int64
    Gens::GAP.GapObj
end

function GapGroupInfo(G::GAP.GapObj)
    return GapGroupInfo(G, Int64(g.Order(G)), g.GeneratorsOfGroup(G))
end

mutable struct LocalCliffordGates
    PauliGates::GapGroupInfo
    NormalizerGates::GapGroupInfo
end

function LocalCliffordGates(n::Int)
    PauliGates = GapGroupInfo(Pauli_group(n))
    NormalizerGates = GapGroupInfo(g.Normalizer(local_clifford_unitaries(n), PauliGates.Grp))
    return LocalCliffordGates(PauliGates, NormalizerGates)
end

function induced_local_clifford_group(CG::LocalCliffordGates)
    P = CG.PauliGates.Grp
    C = CG.NormalizerGates.Grp

    hom = g.ActionHomomorphism(C, P)
    Cl = hom(C)

    return GapGroupInfo(Cl), hom
end

function pauli_to_action_domain_mappings(P::PauliString, G::GAP.GapObj)
    bit_strings = P.bit_strings
    bit_to_int = P.bit_to_int
    group_elements = [f for f in g.Elements(G)]
    n = length(bit_strings[1]) ÷ 2

    internal_labels = Int[]
    action_labels = Int[]

    for a in bit_strings
        for s in 0:1
            signed_pauli = (-1)^s * pauli_n(a)
            internal_label = bit_to_int[a] + s * 2^(2*n)
            action_label = findfirst(w -> w == signed_pauli, group_elements)

            if action_label === nothing
                error("Could not find signed Pauli in GAP group for a = $a, s = $s.")
            end

            push!(internal_labels, internal_label)
            push!(action_labels, action_label)
        end
    end

    return Dict(zip(internal_labels, action_labels)),
           Dict(zip(action_labels, internal_labels))
end

mutable struct LocalCliffordGroup
    n::Int
    PauliInfo::PauliString
    NormalizerGates::GapGroupInfo
    PauliGates::GapGroupInfo
    LocalCliffGroup::GapGroupInfo
    GateGroupHom::GAP.GapObj
    IntActionDict::Dict{Int,Int}
    ActionIntDict::Dict{Int,Int}
end

function LocalCliffordGroup(n::Int)
    CG = LocalCliffordGates(n)
    PauliInfo = PauliString(n)

    LocalCliffGroup, GateGroupHom = induced_local_clifford_group(CG)
    IntActionDict, ActionIntDict = pauli_to_action_domain_mappings(PauliInfo, CG.PauliGates.Grp)

    return LocalCliffordGroup(
        n,
        PauliInfo,
        CG.NormalizerGates,
        CG.PauliGates,
        LocalCliffGroup,
        GateGroupHom,
        IntActionDict,
        ActionIntDict,
    )
end

# -----------------------------------------------------------------------------
# Bloch-vector orbit routines
# A Bloch vector C has length 4^n and represents A = 2^{-n} sum_i C[i] T_i.
# Negative coefficients are encoded by shifting the index by length(C).
# -----------------------------------------------------------------------------

function grouping(C::Vector{<:Number})
    nonzero_indices = findall(!iszero, C)
    nonzero_coefficients = abs.(C[nonzero_indices])

    values = unique(nonzero_coefficients)
    sort!(values)
    value_to_block = Dict(v => i for (i, v) in enumerate(values))

    subsets = [Set{Int}() for _ in values]
    coefficient_dict = Dict{Real,Int}()
    N = length(C)

    for (i, coeff) in enumerate(nonzero_coefficients)
        block = value_to_block[coeff]
        idx = nonzero_indices[i]
        signed_idx = C[idx] >= 0 ? idx : idx + N
        push!(subsets[block], signed_idx)
        coefficient_dict[coeff] = block
    end

    return values, subsets, coefficient_dict
end

function subset_to_vector(N::Int, S::Vector{Vector{Int64}}, dict::Dict{Real,Int})
    block_to_value = reverse_map(dict)
    V = zeros(Real, N)

    for i in 1:length(S)
        abs_val = block_to_value[i]
        for s in S[i]
            if s <= N
                V[s] = abs_val
            else
                V[s - N] = -abs_val
            end
        end
    end

    return V
end

function local_clifford_orbit_of_point(CG::LocalCliffordGroup, C::Vector)
    C = Vector{Real}(C)
    fdict = CG.IntActionDict
    bdict = CG.ActionIntDict
    G = CG.LocalCliffGroup.Grp
    N = length(C)

    _, blocks, value_map = grouping(C)

    S = g.Set(gobj([]))
    for block in blocks
        gap_block = g.Set(gobj([]))
        for e in block
            g.Add(gap_block, fdict[e])
        end
        g.Add(S, g.Set(gap_block))
    end

    S = g.Set(gobj([S]))
    Orbs = g.Orbits(G, S, g.OnSetsSets)[1]

    VOrb = Vector{Vector{Real}}()
    for orb in Orbs
        orbit = Vector{Vector{Int64}}()
        for s in orb
            subset = Int64[bdict[elem] for elem in s]
            push!(orbit, sort(subset))
        end
        push!(VOrb, subset_to_vector(N, orbit, value_map))
    end

    return VOrb
end

function local_clifford_gate_action_on_point(gens, CG::LocalCliffordGroup, C::Vector)
    C = Vector{Real}(C)
    fdict = CG.IntActionDict
    bdict = CG.ActionIntDict
    G = CG.LocalCliffGroup.Grp
    N = length(C)

    group_gens = g.GeneratorsOfGroup(G)
    gg = group_gens[gens[1]]
    for i in gens[2:end]
        gg = gg * group_gens[i]
    end

    _, blocks, value_map = grouping(C)

    image_blocks = Vector{Vector{Int64}}()
    for block in blocks
        image_block = Int64[]
        for e in block
            y = fdict[e]
            gy = g.OnPoints(gobj(y), gg)
            gx = bdict[gy]
            push!(image_block, gx)
        end
        push!(image_blocks, image_block)
    end

    return subset_to_vector(N, image_blocks, value_map)
end

end # module
