using StateVecSim

function check_statevector(circuit, numqubits::Integer)
    statevector = zero_qreg(numqubits)
    apply!(circuit, statevector)
    return statevector
end

@testset "Random Circuit with nearest neighbour gates" begin
    numqubits = 5  # number of qubits for testing
    numcircuits = 4  # number of random circuits to test
    depth = 3
    for test_id in 1:numcircuits
    @debug "Random Circuit #$(test_id)"
        circuit = generate_circuit(RandomCircuitNN(depth), numqubits)
        mps_sv = QbitMPS.compute_statevector(circuit, numqubits)
        check_sv = check_statevector(circuit, numqubits)
        @test real(mps_sv) ≈ check_sv.real && imag(mps_sv) ≈ check_sv.imag
    end
end
