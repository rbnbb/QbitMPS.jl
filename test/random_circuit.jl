@testset "Random Circuit with nearest neighbour gates" begin
    numqubits = 4  # number of qubits for testing
    numcircuits = 2  # number of random circuits to test
    depth = 2
    for test_id in 1:numcircuits
        @debug "Random Circuit #$(test_id)"
        circuit = generate_circuit(RandomCircuitNN(depth), numqubits)
        mps_sv = mps2statevector(simulate_circuit_imps(circuit))
        check_sv = check_statevector(circuit, numqubits)
        @test check_sv ≈ mps_sv
    end
end
