########################################################
# NS(2,3,2) in 15 expectation coordinates
#
# Coordinates:
#   A1,A2,A3, B1,B2,B3, E11,E12,E13,E21,E22,E23,E31,E32,E33
#
# Positivity inequalities:
#   1 + a*A_x + b*B_y + a*b*E_xy >= 0
# for x,y in {1,2,3}, a,b in {+1,-1}.
#
# This gives 36 inequalities in R^15.
#
# Output:
#   - number of vertices
#   - number of combinatorial automorphism orbits
#   - LaTeX table with orbit representatives and orbit sizes
########################################################

using Pkg

base_dir = dirname(@__DIR__)
project_folder = joinpath(base_dir, "MyProject")
Pkg.activate(project_folder)

using Polymake
const pm = Polymake
using GAP
const g = GAP.Globals

########################################################
# Indexing helpers (copied/adapted from symmetries.jl)
########################################################

# Polymake uses 0-based indexing, GAP and Julia use 1-based indexing.
to_one_based_indexing(n::Number) = n + one(n)
to_zero_based_indexing(n::Number) =
    (n > zero(n) ? n - one(n) : throw(ArgumentError("Can't use negative index")))

for f in [:to_one_based_indexing, :to_zero_based_indexing]
    @eval begin
        $f(itr) = $f.(itr)
        $f(s::S) where S<:AbstractSet = S($f.(s))
    end
end

########################################################
# Coordinate bookkeeping
########################################################

const COORD_NAMES = [
    "I",
    "X_1", "Y_1", "Z_1",
    "X_2", "Y_2", "Z_2",
    "XX", "XY", "XZ",
    "YX", "YY", "YZ",
    "ZX", "ZY", "ZZ"
]

# Return the column index in the H-representation row.
# Column 1 is the constant term.
A_col(x::Int) = 1 + x
B_col(y::Int) = 1 + 3 + y
E_col(x::Int, y::Int) = 1 + 6 + 3*(x-1) + y

########################################################
# Build the (2,3,2) nonsignaling polytope
########################################################

function ns232_inequalities()
    rows = Vector{Vector{Int}}()

    # Inequality rows have length 16:
    # [const, A1,A2,A3, B1,B2,B3, E11,...,E33]
    for x in 1:3
        for y in 1:3
            for a in (-1, 1)
                for b in (-1, 1)
                    row = zeros(Int, 16)
                    row[1] = 1
                    row[A_col(x)] = a
                    row[B_col(y)] = b
                    row[E_col(x, y)] = a*b
                    push!(rows, row)
                end
            end
        end
    end

    return reduce(vcat, permutedims.(rows))
end

function ns232_polytope()
    ineq = ns232_inequalities()
    return pm.polytope.Polytope(INEQUALITIES = ineq)
end

########################################################
# Helpers for vertices / representatives / orbits
# (adapted from symmetries.jl)
########################################################

# Input: polymake polytope object P
# Output: GAP group acting on vertices
function vertex_action_GAP(P)
    pm.polytope.combinatorial_symmetries(P)
    pm.give(P, "GROUP")
    gens_poly = P.GROUP.VERTICES_ACTION.GENERATORS
    gens_julia = Vector{Int}.(pm.to_one_based_indexing(gens_poly))
    gens_gap = GAP.Globals.PermList.(GAP.GapObj.(gens_julia))
    return GAP.Globals.Group(gens_gap...)
end

# Input: polymake polytope object
# Output: representative vertices returned by polymake
function representative_vertices(P)
    pm.polytope.combinatorial_symmetries(P)
    pm.give(P, "GROUP.REPRESENTATIVE_VERTICES")
    return P.GROUP.REPRESENTATIVE_VERTICES
end

# Find the row index of vertex v in the vertex matrix V
function find_vertex(v, V)
    nrows = size(V, 1)
    for i in 1:nrows
        if collect(V[i, :]) == collect(v)
            return i
        end
    end
    error("Vertex not found in VERTICES.")
end

# Vertex orbits under the combinatorial automorphism group
function vertex_orbits(P)
    G = vertex_action_GAP(P)
    reps = representative_vertices(P)
    rep_idx = [find_vertex(reps[i, :], P.VERTICES) for i in 1:size(reps, 1)]
    orbit_arr = [sort(Vector{Int}(GAP.Globals.Orbit(G, i))) for i in rep_idx]
    return orbit_arr
end

# Orbit representatives as actual rows of the vertex matrix
function orbit_representatives(P)
    orbits = vertex_orbits(P)
    reps = [collect(P.VERTICES[orb[1], :]) for orb in orbits]
    return reps, orbits
end

########################################################
# Formatting utilities
########################################################

# Convert polymake/GAP rational-ish entries to strings
coord_to_string(x) = string(x)

function vector_to_latex(v; names=COORD_NAMES)
    parts = String[]
    for i in 1:length(v)
        push!(parts, "$(names[i]) = $(coord_to_string(v[i]))")
    end
    return "\\(" * join(parts, ",\\; ") * "\\)"
end

function short_vector_to_latex(v)
    return "\\((" * join(coord_to_string.(v), ", ") * ")\\)"
end

########################################################
# LaTeX output
########################################################

function write_orbit_table(filename::AbstractString, reps, orbits)
    open(filename, "w") do io
        println(io, "\\begin{table}[ht]")
        println(io, "\\centering")
        println(io, "\\begin{tabular}{c c p{0.72\\linewidth}}")
        println(io, "\\hline")
        println(io, "Orbit & Size & Representative in expectation coordinates\\\\")
        println(io, "\\hline")

        for (k, rep) in enumerate(reps)
            rep_tex = short_vector_to_latex(rep)
            println(io, "\$\\mathcal{O}_{$k}\$ & $(length(orbits[k])) & $rep_tex\\\\")
        end

        println(io, "\\hline")
        println(io, "\\end{tabular}")
        println(io, "\\caption{Orbits of the vertices of the \$(2,3,2)\$ nonsignaling polytope in the 15 expectation coordinates under the action of the combinatorial automorphism group.}")
        println(io, "\\label{tab:ns232_orbits}")
        println(io, "\\end{table}")
    end
end

function write_verbose_orbit_table(filename::AbstractString, reps, orbits)
    open(filename, "w") do io
        println(io, "\\begin{table}[ht]")
        println(io, "\\centering")
        println(io, "\\footnotesize")
        println(io, "\\setlength{\\tabcolsep}{6pt}")
        println(io, "\\renewcommand{\\arraystretch}{1.2}")
        println(io, "\\begin{tabular}{c c cccccccccccccccc}")
        println(io, "\\hline")

        # Header row
        header = ["Orbit", "Size", "\$I\$"]
        append!(header, ["\$X_1\$", "\$Y_1\$", "\$Z_1\$",
                         "\$X_2\$", "\$Y_2\$", "\$Z_2\$",
                         "\$XX\$", "\$XY\$", "\$XZ\$",
                         "\$YX\$", "\$YY\$", "\$YZ\$",
                         "\$ZX\$", "\$ZY\$", "\$ZZ\$"])
        println(io, join(header, " & ") * " \\\\")
        println(io, "\\hline")

        # Data rows
        for (k, rep) in enumerate(reps)
            vals = string.(rep)
            row = ["\$\\mathcal{O}_{$k}\$", string(length(orbits[k]))]
            append!(row, vals)
            println(io, join(row, " & ") * " \\\\")
        end

        println(io, "\\hline")
        println(io, "\\end{tabular}")
        println(io, "\\caption{Orbit representatives for the vertices of the \$(2,3,2)\$ nonsignaling polytope in expectation coordinates.}")
        println(io, "\\label{tab:ns232_orbits_verbose}")
        println(io, "\\end{table}")
    end
end

########################################################
# Main
########################################################

function main()
    println("Building the (2,3,2) nonsignaling polytope in 15 expectation coordinates...")
    P = ns232_polytope()

    # Force computation of vertices if needed
    pm.give(P, "VERTICES")

    nv = size(P.VERTICES, 1)
    println("Number of vertices: ", nv)

    reps, orbits = orbit_representatives(P)

    println("Number of combinatorial automorphism orbits: ", length(orbits))
    println("Orbit sizes: ", [length(O) for O in orbits])

    # Expect 5 according to the manuscript discussion
    if length(orbits) != 5
        @warn "Expected 5 orbits, but found $(length(orbits))."
    end

    write_orbit_table("ns232_orbits_table.tex", reps, orbits)
    write_verbose_orbit_table("ns232_orbits_table_verbose.tex", reps, orbits)

    println("Wrote:")
    println("  ns232_orbits_table.tex")
    println("  ns232_orbits_table_verbose.tex")

    println("\nOrbit representatives:")
    for (k, rep) in enumerate(reps)
        println("Orbit ", k, " (size ", length(orbits[k]), "): ", rep)
    end
end

main()