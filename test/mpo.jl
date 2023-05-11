using QbitMPS: simulate_circuit_mpo, random_circuit_lr, random_circuit_nn, mps2statevector
using QbitMPS: Circuit, GateH, GateCX, GateU

@testset "Test 1 qubit gates for MPO simulator" begin
    function test_1qubit_unitaries(numqubits::Integer)
        function random_circuit_only_1q_unitaries(numqubits::Integer)
            circuit = Circuit()
            for qubit in 1:numqubits
                push!(circuit, GateU(rand(3)...), qubit)
            end
            return circuit
        end
        @assert numqubits > 1
        circ = random_circuit_only_1q_unitaries(numqubits)
        mps_sv = mps2statevector(simulate_circuit_mpo(circ))
        sv = check_statevector(circ, max(4, numqubits))[1:2^numqubits]
        # stringify(ind) = [string(string.(j...)) for j in ind]
        # return [stringify(list_fock_states(numqubits)) round.(sv; digits = 3) round.( mps_sv; digits = 3,)]
        return mps_sv ≈ sv
    end
    MAXQUBITS = 7
    for n in 3:MAXQUBITS
        @test test_1qubit_unitaries.(n)
    end
end

@testset "Test 2 qubit nearest neighbour gates for MPO simulator" begin
    function test_bell_preparation_with_CX()
        circ = Circuit()
        push!(circ, GateH(), 1)
        push!(circ, GateCX(), 1, 2)
        return mps2statevector(simulate_circuit_mpo(circ)) ≈ check_statevector(circ, 4)[1:4]
    end
    function test_2q_nn(numqubits::Integer, depth::Integer)
        @assert numqubits > 1
        circ = random_circuit_nn(numqubits, depth)
        mps_sv = mps2statevector(simulate_circuit_mpo(circ))
        sv = check_statevector(circ, max(4, numqubits))[1:2^numqubits]
        return mps_sv ≈ sv
    end
    @test test_bell_preparation_with_CX()
    MAXQUBITS = 5
    CIRCUIT_DEPTH = 2
    for n in 2:MAXQUBITS
        @test test_2q_nn(n, CIRCUIT_DEPTH)
    end
end

@testset "Test 2 qubits long-range gates for MPO simulator" begin
    function test_bell_preparation_3q()
        circ = Circuit()
        push!(circ, GateH(), 1)
        push!(circ, GateCX(), 1, 3)
        mps_sv = mps2statevector(simulate_circuit_mpo(circ))
        sv = check_statevector(circ, 4)[1:2^3]
        # stringify(ind) = [string(string.(j...)) for j in ind]
        # return [stringify(list_fock_states(3)) round.(sv; digits = 3) round.( mps_sv; digits = 3,)]
        return mps_sv ≈ sv
    end
    function test_2q_long_range(numqubits::Integer, depth::Integer)
        @assert numqubits > 1
        circ = random_circuit_lr(numqubits, depth)
        mps_sv = mps2statevector(simulate_circuit_mpo(circ))
        sv = check_statevector(circ, max(4, numqubits))[1:2^numqubits]
        return mps_sv ≈ sv
    end
    @test test_bell_preparation_3q()
    MAXQUBITS = 5
    CIRCUIT_DEPTH = 2
    for n in 2:MAXQUBITS
        @test test_2q_long_range(n, CIRCUIT_DEPTH)
    end
end
