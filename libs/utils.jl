
"""
    all_dit_strings(d, N) -> Vector{Vector{Int}}

Generate all possible dit strings of length `N`, where each element is an integer in the range `0:d-1`.

# Parameters
- `d::Int`: The number of outcomes (each digit will be in the range 0 to d-1).
- `N::Int`: The length of each dit string.

# Returns
A vector containing all dit strings, each represented as a vector of integers.
"""
function all_dit_strings(d, N)
    return [collect(t) for t in Iterators.product(ntuple(_ -> 0:d-1, N)...)]
end


function array_to_matrix(A)
    n = length(A[1]); N = length(A);
    R = Array{Rational{Int64}}(undef,n,0)
    for i in 1:N
        R = hcat(R,A[i]);
    end
    return R
end

# example usage: save_input_lrs(A2,"L2.ine","Lambda2","Lambda2"," ","H")