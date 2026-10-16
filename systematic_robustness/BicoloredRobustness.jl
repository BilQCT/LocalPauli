module BicoloredRobustness

using GLPK
using HDF5
using JSON3
using JuMP
using LinearAlgebra
using SHA
include(joinpath(@__DIR__, "..", "libs", "MaxWeightSeed.jl"))

const MOI = JuMP.MOI

export BicoloredRecord,
       FixedDictionary,
       RobustnessResult,
       BatchResult,
       bicolored_state,
       generate_bicolored_states,
       load_bicolored_records,
       load_fixed_dictionary,
       compute_robustness,
       compute_batch,
       save_batch,
       validate_pauli_state

"""One color-preserving graph-isomorphism representative from the catalogue."""
struct BicoloredRecord
    id::Int
    n::Int
    edges::Vector{Tuple{Int,Int}}
    red::Vector{Int}
    hash::String
end

"""A fixed collection of input vectors, stored as columns in Pauli coordinates."""
struct FixedDictionary{T<:Real}
    vectors::Matrix{T}
    n::Int
    label::String
    sources::Vector{String}
end

"""Numerical result of one fixed-dictionary robustness linear program."""
struct RobustnessResult
    objective::Float64
    primal_residual_inf::Float64
    trace_residual::Float64
    termination_status::String
    primal_status::String
    support_indices::Vector{Int}
    support_coefficients::Vector{Float64}
    dual::Vector{Float64}
    dual_violation::Float64
    duality_gap::Float64
end

"""Results for a collection of targets evaluated against one fixed dictionary."""
struct BatchResult
    dictionary_label::String
    dictionary_sources::Vector{String}
    dictionary_columns::Int
    n::Int
    records::Vector{BicoloredRecord}
    results::Vector{RobustnessResult}
end

const LOCAL_I = 0
const LOCAL_X = 1
const LOCAL_Y = 2
const LOCAL_Z = 3

"""Convert lexicographic local Pauli labels (I,X,Y,Z = 0,1,2,3) to an index."""
function pauli_index(labels::AbstractVector{<:Integer})
    index0 = 0
    for label in labels
        0 <= label <= 3 || throw(ArgumentError("local Pauli label must lie in 0:3"))
        index0 = 4 * index0 + label
    end
    return index0 + 1
end

"""Return the lexicographic local Pauli labels for a one-based vector index."""
function pauli_labels(index::Int, n::Int)
    1 <= index <= 4^n || throw(BoundsError(1:4^n, index))
    value = index - 1
    labels = zeros(Int, n)
    for q in n:-1:1
        labels[q] = value % 4
        value ÷= 4
    end
    return labels
end

local_xz(label::Int) = label == LOCAL_I ? (0, 0) :
                       label == LOCAL_X ? (1, 0) :
                       label == LOCAL_Y ? (1, 1) : (0, 1)

function xz_label(x::Int, z::Int)
    x == 0 && z == 0 && return LOCAL_I
    x == 1 && z == 0 && return LOCAL_X
    x == 1 && z == 1 && return LOCAL_Y
    return LOCAL_Z
end

"""Apply a signed permutation of the Pauli basis to a coefficient vector."""
function apply_signed_permutation(map_label, coefficients::AbstractVector{<:Real})
    output = zeros(Float64, length(coefficients))
    n_float = log(4, length(coefficients))
    n = round(Int, n_float)
    4^n == length(coefficients) || throw(DimensionMismatch("Pauli vector length is not 4^n"))

    for index in eachindex(coefficients)
        coefficient = coefficients[index]
        iszero(coefficient) && continue
        labels, sign = map_label(pauli_labels(index, n))
        output[pauli_index(labels)] += sign * coefficient
    end
    return output
end

function apply_hadamard(coefficients::AbstractVector{<:Real}, q::Int)
    return apply_signed_permutation(coefficients) do labels
        mapped = copy(labels)
        label = labels[q]
        sign = label == LOCAL_Y ? -1.0 : 1.0
        mapped[q] = label == LOCAL_X ? LOCAL_Z : label == LOCAL_Z ? LOCAL_X : label
        return mapped, sign
    end
end

function apply_cz(coefficients::AbstractVector{<:Real}, control::Int, target::Int)
    control == target && throw(ArgumentError("CZ endpoints must be distinct"))
    return apply_signed_permutation(coefficients) do labels
        mapped = copy(labels)
        xc, zc = local_xz(labels[control])
        xt, zt = local_xz(labels[target])
        mapped[control] = xz_label(xc, (zc + xt) % 2)
        mapped[target] = xz_label(xt, (zt + xc) % 2)
        phase = (xc * xt * zc + xc * xt * zt) % 2
        return mapped, phase == 0 ? 1.0 : -1.0
    end
end

function apply_t(coefficients::AbstractVector{<:Real}, q::Int)
    output = zeros(Float64, length(coefficients))
    n = round(Int, log(4, length(coefficients)))
    invsqrt2 = inv(sqrt(2.0))

    for index in eachindex(coefficients)
        coefficient = coefficients[index]
        iszero(coefficient) && continue
        labels = pauli_labels(index, n)
        label = labels[q]
        if label == LOCAL_I || label == LOCAL_Z
            output[index] += coefficient
        elseif label == LOCAL_X
            output[index] += coefficient * invsqrt2
            labels[q] = LOCAL_Y
            output[pauli_index(labels)] += coefficient * invsqrt2
        else
            output[index] += coefficient * invsqrt2
            labels[q] = LOCAL_X
            output[pauli_index(labels)] -= coefficient * invsqrt2
        end
    end
    return output
end

"""
    bicolored_state(n, edges, red)

Construct the Pauli vector `Tr(P rho)` of

    |psi_G,U> = T^⊗n H^⊗U E(G) |+>^⊗n.

Qubits and edge endpoints are one-based.
"""
function bicolored_state(n::Int,
                         edges::AbstractVector{<:Tuple{Int,Int}},
                         red::AbstractVector{<:Integer})
    n >= 1 || throw(ArgumentError("n must be positive"))
    all(q -> 1 <= q <= n, red) || throw(ArgumentError("red vertex outside 1:n"))
    length(unique(red)) == length(red) || throw(ArgumentError("red vertices must be unique"))

    coefficients = [1.0]
    for _ in 1:n
        coefficients = kron(coefficients, [1.0, 1.0, 0.0, 0.0])
    end

    for (i, j) in edges
        1 <= i <= n && 1 <= j <= n || throw(ArgumentError("edge endpoint outside 1:n"))
        coefficients = apply_cz(coefficients, i, j)
    end
    for q in red
        coefficients = apply_hadamard(coefficients, q)
    end
    for q in 1:n
        coefficients = apply_t(coefficients, q)
    end

    validate_pauli_state(coefficients; pure=true)
    return coefficients
end

bicolored_state(record::BicoloredRecord) =
    bicolored_state(record.n, record.edges, record.red)

"""Check normalization, finiteness, and optionally the pure-state Pauli norm."""
function validate_pauli_state(coefficients::AbstractVector{<:Real};
                              pure::Bool=false,
                              atol::Float64=1e-9)
    all(isfinite, coefficients) || throw(ArgumentError("Pauli vector contains non-finite values"))
    n = round(Int, log(4, length(coefficients)))
    4^n == length(coefficients) || throw(DimensionMismatch("Pauli vector length is not 4^n"))
    isapprox(coefficients[1], 1.0; atol=atol, rtol=0) ||
        throw(ArgumentError("identity coefficient must equal one"))
    if pure
        isapprox(sum(abs2, coefficients), 2.0^n; atol=atol * 2.0^n, rtol=0) ||
            throw(ArgumentError("Pauli vector does not have the norm of a pure state"))
    end
    return true
end

"""Load the zero-based JSON catalogue and return one-based Julia records."""
function load_bicolored_records(path::AbstractString)
    raw_records = JSON3.read(read(path, String))
    records = BicoloredRecord[]
    for raw in raw_records
        edges = Tuple{Int,Int}[(Int(e[1]) + 1, Int(e[2]) + 1) for e in raw["edges"]]
        red = Int[Int(q) + 1 for q in raw["red"]]
        hash = haskey(raw, "hash") ? String(raw["hash"]) : ""
        push!(records, BicoloredRecord(Int(raw["id"]), Int(raw["n"]), edges, red, hash))
    end
    return records
end

"""Generate a matrix whose columns follow the catalogue order."""
function generate_bicolored_states(records::AbstractVector{BicoloredRecord})
    isempty(records) && throw(ArgumentError("record collection is empty"))
    n = only(unique(record.n for record in records))
    states = Matrix{Float64}(undef, 4^n, length(records))
    for (column, record) in enumerate(records)
        states[:, column] = bicolored_state(record)
    end
    return states
end

function read_hdf5_matrix(path::AbstractString, dataset::AbstractString)
    MaxWeightSeed.validate_maxweight_file(path)
    return h5open(path, "r") do file
        haskey(file, dataset) || throw(KeyError("dataset '$dataset' not found in $path"))
        read(file[dataset])
    end
end

"""
    load_fixed_dictionary(paths; dataset="data", label, deduplicate=false)

Load one or more fixed input-vector matrices and concatenate their columns.
Every matrix must use the lexicographic I,X,Y,Z Pauli convention and have
identity row equal to one.
"""
function load_fixed_dictionary(paths::AbstractVector{<:AbstractString};
                               dataset::AbstractString="data",
                               label::AbstractString=join(basename.(paths), "+"),
                               deduplicate::Bool=false)
    isempty(paths) && throw(ArgumentError("at least one dictionary path is required"))
    blocks = [read_hdf5_matrix(path, dataset) for path in paths]
    row_count = only(unique(size(block, 1) for block in blocks))
    n = round(Int, log(4, row_count))
    4^n == row_count || throw(DimensionMismatch("dictionary row count is not 4^n"))
    all(block -> size(block, 1) == row_count, blocks) ||
        throw(DimensionMismatch("dictionary blocks have different row counts"))

    vectors = Matrix{Float64}(hcat(blocks...))
    all(isfinite, vectors) || throw(ArgumentError("dictionary contains non-finite entries"))
    all(isapprox.(vectors[1, :], 1.0; atol=1e-12, rtol=0)) ||
        throw(ArgumentError("dictionary identity row must equal one"))

    if deduplicate
        seen = Set{Vector{Float64}}()
        keep = Int[]
        for column in axes(vectors, 2)
            value = vectors[:, column]
            if !(value in seen)
                push!(seen, copy(value))
                push!(keep, column)
            end
        end
        vectors = vectors[:, keep]
    end
    return FixedDictionary(vectors, n, String(label), String.(paths))
end

load_fixed_dictionary(path::AbstractString; kwargs...) =
    load_fixed_dictionary([path]; kwargs...)

"""Compute robustness relative to exactly the columns in `dictionary`."""
function compute_robustness(dictionary::FixedDictionary,
                            target::AbstractVector{<:Real};
                            support_tolerance::Float64=1e-9,
                            residual_tolerance::Float64=1e-7,
                            silent::Bool=true)
    length(target) == size(dictionary.vectors, 1) ||
        throw(DimensionMismatch("target and dictionary dimensions differ"))
    validate_pauli_state(target)
    all(isfinite, dictionary.vectors) || error("Non-finite dictionary")
    all(isapprox.(dictionary.vectors[1, :], 1; atol=1e-12, rtol=0)) ||
        error("Dictionary identity coordinates must equal one")

    model = Model(GLPK.Optimizer)
    silent && set_silent(model)
    number_of_vectors = size(dictionary.vectors, 2)
    @variable(model, coefficients[1:number_of_vectors])
    @variable(model, magnitudes[1:number_of_vectors] >= 0)
    @objective(model, Min, sum(magnitudes))
    @constraint(model, magnitudes .>= coefficients)
    @constraint(model, magnitudes .>= -coefficients)
    reconstruction = @constraint(model, dictionary.vectors * coefficients .== target)
    optimize!(model)

    termination = termination_status(model)
    termination == MOI.OPTIMAL || error("robustness LP terminated with $termination")
    values = value.(coefficients)
    residual = dictionary.vectors * values - target
    residual_inf = norm(residual, Inf)
    trace_residual = abs(sum(values) - target[1])
    residual_inf <= residual_tolerance ||
        error("LP reconstruction residual $residual_inf exceeds tolerance $residual_tolerance")
    support = findall(value -> abs(value) > support_tolerance, values)
    # Validate the actual saved sparse decomposition, not just the untrimmed LP.
    saved_residual = norm(dictionary.vectors[:, support] * values[support] - target, Inf)
    saved_residual <= residual_tolerance || error("Truncated decomposition fails reconstruction")
    witness = dual.(reconstruction)
    dual_violation = max(0.0, maximum(abs.(dictionary.vectors' * witness)) - 1)
    gap = abs(sum(abs, values[support]) - dot(target, witness))
    max(dual_violation, gap) <= residual_tolerance || error("Primal/dual certificate failed")

    return RobustnessResult(
        objective_value(model),
        max(residual_inf, saved_residual),
        trace_residual,
        string(termination),
        string(primal_status(model)),
        support,
        values[support],
        witness,
        dual_violation,
        gap,
    )
end

"""Evaluate catalogue records against one dictionary, without changing it."""
function compute_batch(dictionary::FixedDictionary,
                       records::AbstractVector{BicoloredRecord};
                       kwargs...)
    isempty(records) && throw(ArgumentError("record collection is empty"))
    all(record -> record.n == dictionary.n, records) ||
        throw(DimensionMismatch("record qubit count and dictionary qubit count differ"))
    results = [compute_robustness(dictionary, bicolored_state(record); kwargs...)
               for record in records]
    return BatchResult(dictionary.label, dictionary.sources,
                       size(dictionary.vectors, 2), dictionary.n,
                       collect(records), results)
end

"""Save objectives, diagnostics, catalogue metadata, and sparse decompositions."""
function save_batch(path::AbstractString, batch::BatchResult)
    mkpath(dirname(abspath(path)))
    h5open(path, "w") do file
        attributes(file)["schema_version"] = 2
        attributes(file)["state_family"] = "T^n H^U E(G) |+>^n"
        attributes(file)["pauli_convention"] = "lexicographic I,X,Y,Z; entries Tr(P rho)"
        attributes(file)["dictionary_label"] = batch.dictionary_label
        attributes(file)["dictionary_columns"] = batch.dictionary_columns
        attributes(file)["n"] = batch.n
        file["dictionary_sources"] = batch.dictionary_sources
        file["dictionary_sha256"] = [bytes2hex(open(sha256, source)) for source in batch.dictionary_sources]
        file["record_id"] = [record.id for record in batch.records]
        file["record_hash"] = [record.hash for record in batch.records]
        file["red_vertices"] = [join(record.red, ",") for record in batch.records]
        file["edges"] = [join(("$(i)-$(j)" for (i, j) in record.edges), ",")
                         for record in batch.records]
        file["objective"] = [result.objective for result in batch.results]
        file["primal_residual_inf"] = [result.primal_residual_inf for result in batch.results]
        file["trace_residual"] = [result.trace_residual for result in batch.results]
        file["dual_violation"] = [result.dual_violation for result in batch.results]
        file["duality_gap"] = [result.duality_gap for result in batch.results]
        file["termination_status"] = [result.termination_status for result in batch.results]
        file["primal_status"] = [result.primal_status for result in batch.results]
        decompositions = create_group(file, "decompositions")
        for (index, result) in enumerate(batch.results)
            group = create_group(decompositions, string(index))
            group["indices"] = result.support_indices
            group["coefficients"] = result.support_coefficients
            group["dual"] = result.dual
        end
    end
    return path
end

end # module
