using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "local_robustness", "robustness_functions.jl"))

@testset "legacy helpers use physical expectation coordinates" begin
    @test pauli([1,1]) == ComplexF64[0 -im; im 0]
    @test norm(plus_state(3)) ≈ 1
    @test convert_matrix_to_pauli_basis(1, ComplexF64[1 0; 0 0]) ≈ [1,0,0,1]
    for n in 2:3
        ps = PauliString(n)
        states, _ = generate_graph_state_vectors(n)
        for bits in all_dit_strings(2,n)
            ket=states[1]
            for q in 1:n
                if bits[q]==1
                    ket=hadamard_gate(n,q)*ket
                end
                ket=t_gate(n,q)*ket
            end
            expected=BicoloredRobustness.bicolored_state(n,[(i,i+1) for i in 1:n-1],findall(==(1),bits))
            @test convert_matrix_to_pauli_basis(n,ket*ket') ≈ expected atol=1e-12
        end
    end
end
