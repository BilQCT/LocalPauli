function reverse_map(map)
    in_key = collect(keys(map)); out_value = collect(values(map))
    return Dict(zip(out_value,in_key))
end


function composition(map1,map2)
    Yin = Set(values(map1)); Yout = Set(keys(map2))

    if isempty(intersect!(Yin,Yout)) == true
        error("Codomain of first map have elements in domain of second.")
    else
        Y = intersect!(Set(values(map1)), Set(keys(map2)))
        in_key = [reverse_map(map1)[y] for y in Y]
        out_value = [map2[y] for y in Y]
        return Dict(zip(in_key,out_value))
    end
end



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