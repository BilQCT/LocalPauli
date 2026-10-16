using Pkg

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(joinpath(REPOSITORY_ROOT, "MyProject"))

include(joinpath(@__DIR__, "BicoloredRobustness.jl"))
using .BicoloredRobustness

function usage()
    println("""
    Usage:
      julia --project=MyProject systematic_robustness/run_batch.jl \\
        N OUTPUT.h5 LABEL INPUT_1.h5 [INPUT_2.h5 ...]

    Example (a nested four-qubit comparison dictionary):
      julia --project=MyProject systematic_robustness/run_batch.jl \\
        4 results/L2_n4.h5 L2 \\
        data/stabilizer_states_4.h5 data/deterministic_4.h5 \\
        data/max_weight_operators_p_2_4.h5

    All relative paths are resolved from the LocalPauli repository root.
    """)
end

function repository_path(path::AbstractString)
    return isabspath(path) ? normpath(path) : normpath(joinpath(REPOSITORY_ROOT, path))
end

function main(arguments)
    length(arguments) >= 4 || begin
        usage()
        return 2
    end

    n = parse(Int, arguments[1])
    output_path = repository_path(arguments[2])
    label = arguments[3]
    dictionary_paths = repository_path.(arguments[4:end])
    catalogue_path = joinpath(REPOSITORY_ROOT, "data", "bicolored_graphs_n$(n).json")

    records = load_bicolored_records(catalogue_path)
    dictionary = load_fixed_dictionary(dictionary_paths; label=label, deduplicate=false)

    println("State family: T^⊗n H^⊗U E(G) |+>^⊗n")
    println("Catalogue: $catalogue_path ($(length(records)) records)")
    println("Dictionary: $(dictionary.label) ($(size(dictionary.vectors, 2)) fixed vectors)")
    println("Output: $output_path")

    batch = compute_batch(dictionary, records)
    save_batch(output_path, batch)

    objectives = [result.objective for result in batch.results]
    residuals = [result.primal_residual_inf for result in batch.results]
    println("Completed $(length(batch.results)) LPs")
    println("Objective range: $(extrema(objectives))")
    println("Maximum reconstruction residual: $(maximum(residuals))")
    return 0
end

exit(main(ARGS))
