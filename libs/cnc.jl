using LinearAlgebra

function do_commute(T_a::Vector{Int}, T_b::Vector{Int})
    n = length(T_a) ÷ 2
    T_ax = T_a[1:n]
    T_az = T_a[n+1:end]

    T_bx = T_b[1:n]
    T_bz = T_b[n+1:end]

    return (dot(T_az, T_bx) + dot(T_ax, T_bz)) % 2 == 0
end

function find_isotropic_gens(isotropic::Set{Vector{Int}})
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
        newly_gen_list = [a + b for a in generated_set]
        gen_list = append!(newly_gen_list, [a for a in generated_set])
        newly_generated_set = Set{Vector{Int}}(newly_gen_list)
        generated_set = Set{Vector{Int}}(gen_list)
    end

    return isotropic_gens
end

mutable struct CncSet
    isotropic_gens::Set{Vector{Int}}
    anticommuting_paulis::Set{Vector{Int}}
    n::Int
    m::Int

    function CncSet(isotropic_gens::Set{Vector{Int}}, anticommuting_paulis::Set{Vector{Int}})
        n = length(first(isotropic_gens)) ÷ 2
        m = (length(anticommuting_paulis) - 1) ÷ 2
        new(isotropic_gens, anticommuting_paulis, n, m)
    end

    function CncSet(cnc_set::Set{Vector{Int}})
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
            anticommuting_paulis = Set{Vector{Int}}([[0 for i in 1:2n]])
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
            isotropic_gens = Set{Vector{Int}}([[0 for i in 1:2n]])
        end
        
        m = (length(anticommuting_paulis) - 1) ÷ 2
        new(isotropic_gens, anticommuting_paulis, n, m)
    end
end


function show(io::IO, c::CncSet)
    d = Dict()
    d["isotropic_gens"] = c.isotropic_gens
    d["anticommuting_paulis"] = c.anticommuting_paulis
    show(io, d)
end