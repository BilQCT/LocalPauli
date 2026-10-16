using Pkg
using Test
using LinearAlgebra
using HDF5

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(REPOSITORY_ROOT, "MyProject"))

include(joinpath(REPOSITORY_ROOT, "systematic_robustness", "BicoloredRobustness.jl"))
using .BicoloredRobustness

@testset "bicolored state generation" begin
    plus = bicolored_state(1, Tuple{Int,Int}[], Int[])
    @test plus ≈ [1.0, inv(sqrt(2.0)), inv(sqrt(2.0)), 0.0]

    decorated = bicolored_state(1, Tuple{Int,Int}[], [1])
    @test decorated ≈ [1.0, 0.0, 0.0, 1.0]

    bell_graph = bicolored_state(2, [(1, 2)], Int[])
    @test bell_graph[1] ≈ 1.0
    @test sum(abs2, bell_graph) ≈ 4.0
end

@testset "catalogue and stored-state agreement" begin
    for n in 2:4
        records = load_bicolored_records(joinpath(REPOSITORY_ROOT, "data", "bicolored_graphs_n$(n).json"))
        generated = generate_bicolored_states(records)
        stored = h5open(joinpath(REPOSITORY_ROOT, "data", "bicolored_resource_projectors_pauli_n$(n).h5"), "r") do file
            read(file["data"])
        end
        @test generated ≈ stored atol=1e-12
    end
end

@testset "fixed-dictionary robustness" begin
    path = joinpath(REPOSITORY_ROOT, "data", "deterministic_2.h5")
    dictionary = load_fixed_dictionary(path; label="D2")
    records = load_bicolored_records(joinpath(REPOSITORY_ROOT, "data", "bicolored_graphs_n2.json"))
    result = compute_robustness(dictionary, bicolored_state(records[1]))
    @test result.objective ≈ 1.2071067811865475 atol=1e-8
    @test result.primal_residual_inf <= 1e-7
    @test result.trace_residual <= 1e-7
    @test result.dual_violation <= 1e-8
    @test result.duality_gap <= 1e-8
    @test dot(bicolored_state(records[1]), result.dual) ≈ result.objective atol=1e-8
    batch = BatchResult(dictionary.label, dictionary.sources,
                        size(dictionary.vectors, 2), dictionary.n,
                        records[1:1], [result])
    mktempdir() do directory
        output = save_batch(joinpath(directory, "certificate.h5"), batch)
        h5open(output, "r") do file
            @test read(attributes(file)["schema_version"]) == 2
            @test read(file["objective"]) == [result.objective]
            @test read(file["decompositions/1/dual"]) == result.dual
            @test read(file["dictionary_sha256"]) == [bytes2hex(open(BicoloredRobustness.sha256, path))]
            indices = read(file["decompositions/1/indices"])
            coefficients = read(file["decompositions/1/coefficients"])
            @test norm(dictionary.vectors[:, indices] * coefficients - bicolored_state(records[1]), Inf) < 1e-8
        end
    end
end

@testset "manuscript seed arithmetic and stale-file rejection" begin
    seed = BicoloredRobustness.MaxWeightSeed
    @test seed.gamma_elementary_p([1, 1, 1], 2) == 1
    @test seed.gamma_elementary_p([2, 2], 2) == 1 # distinguishes product mod 3 from product mod 2
    @test seed.gamma_elementary_p([1, 1, 1, 0], 2) == 1
    @test seed.gamma_elementary_p([2, 2, 2], 3) == 0
    @test seed.gamma_elementary_p([2, 2, 1], 3) == 1
    @test seed.gamma_elementary_p([1, 1], 3) == 0
    @test_throws ArgumentError seed.gamma_elementary_p([3, 1], 2)
    for n in 2:4
        path = joinpath(REPOSITORY_ROOT, "data", "max_weight_operators_p_2_$(n).h5")
        @test seed.validate_maxweight_file(path)
        h5open(path, "r") do f
            @test read(f["fiducial"]) == seed.seed_coordinates(n)
            @test size(f["data"]) == (4^n, Dict(2=>96, 3=>13824, 4=>165888)[n])
        end
    end
    mktempdir() do directory
        stale = joinpath(directory, "max_weight_operators_p_2_2.h5")
        h5open(stale, "w") do file
            file["data"] = ones(16, 1)
        end
        @test_throws ErrorException load_fixed_dictionary(stale; label="stale")
    end
end

@testset "trace-one CNC dictionary" begin
    d = load_fixed_dictionary(joinpath(REPOSITORY_ROOT, "data", "cnc_operators_3.h5"); label="C")
    @test all(d.vectors[1, :] .== 1)
    @test size(d.vectors) == (64, 71136)
end
