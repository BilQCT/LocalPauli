module MaxWeightSeed

using Combinatorics
using HDF5

export gamma_elementary_p, seed_coordinates, validate_maxweight_file,
       SEED_CONVENTION, ORBIT_CONVENTION

const SEED_CONVENTION = "termwise_product_mod3_sum_mod2_v1"
const ORBIT_CONVENTION = "proper_local_clifford"

"""Manuscript seed: multiply each monomial in Z3, then add its parity in Z2.

There is deliberately NO reduction of the running sum modulo 3.
For example gamma_elementary_p([1,1,1], 2) == 1.
"""
function gamma_elementary_p(u::AbstractVector{<:Integer}, p::Integer)
    p >= 1 || throw(ArgumentError("p must be positive"))
    all(x -> 0 <= x <= 2, u) || throw(ArgumentError("inputs must lie in 0:2"))
    total_z2 = 0
    for indices in combinations(eachindex(u), p)
        product_z3 = 1
        for i in indices
            product_z3 = mod(product_z3 * u[i], 3)
        end
        total_z2 = mod(total_z2 + product_z3, 2)
    end
    return total_z2
end

"""Coordinates Tr(P A), lexicographic I,X,Y,Z, including identity 1."""
function seed_coordinates(n::Integer, p::Integer=2)
    2 <= p <= n || throw(ArgumentError("need 2 <= p <= n"))
    result = zeros(Int8, 4^n)
    result[1] = 1
    for input in Iterators.product(ntuple(_ -> 0:2, n)...)
        index0 = foldl((index, digit) -> 4 * index + digit + 1, input; init=0)
        result[index0 + 1] = iszero(gamma_elementary_p(collect(input), p)) ? 1 : -1
    end
    return result
end

"""Reject old max-weight dictionaries before they can silently enter an LP.

Full membership and completeness are checked by the numerical audit; this
loader guard checks the versioned convention and actual stored fiducial.
"""
function validate_maxweight_file(path::AbstractString)
    occursin("max_weight_operators", basename(path)) || return true
    h5open(path, "r") do file
        attrs = attributes(file)
        for (key, expected) in (("seed_convention", SEED_CONVENTION),
                                ("orbit_convention", ORBIT_CONVENTION))
            haskey(attrs, key) && read(attrs[key]) == expected ||
                error("Stale or unsupported max-weight dictionary: $path. Regenerate it with the manuscript seed.")
        end
        n, p = Int(read(file["n"])), Int(read(file["p"]))
        read(file["fiducial"]) == seed_coordinates(n, p) ||
            error("Stored max-weight fiducial disagrees with the manuscript: $path")
        size(file["data"], 1) == 4^n || error("Incorrect Pauli dimension in $path")
    end
    return true
end

end
