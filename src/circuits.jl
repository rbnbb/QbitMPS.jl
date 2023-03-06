abstract type QuantumAlgorithm end
struct QFT <: QuantumAlgorithm end
struct RandomCircuitNN <: QuantumAlgorithm
    depth::Integer
end

generate_circuit(::Type{QFT}, numqubits::Integer) = quantum_fourier_circuit(numqubits)
generate_circuit(algorithm::RandomCircuitNN, numqubits::Integer) = random_circuit_nn(numqubits, algorithm.depth)

@doc raw"""
    quantum_fourier_circuit(numqubits::Int)

Return the Quantum Fourier Transform circuit for numqubits.

The formula for QFT used is:
using ITensors: random_unitary
using ITensors: randomU
```math
QFT_N | x \rangle = \frac{1}{\sqrt{N}} \sum_{y=0}^{N-1} e^{\frac{2 \pi i xy }{ 2^n}} | y \rangle
```
with ``|x\rangle`` a number state and ``N=2^{\text{numqubits}}``.
"""
function quantum_fourier_circuit(numqubits::Integer)::Circuit
    circuit = Circuit(numqubits)
    for qubit in 1:numqubits
        push!(circuit, GateH(), qubit)
        for target_qubit in qubit+1:numqubits
            phase = 2.0 * π / (2.0^(target_qubit - qubit + 1.0))
            push!(circuit, GateCP(phase), qubit, target_qubit)
        end
    end
    for qubit in 1:numqubits/2
        push!(circuit, GateSWAP(), qubit, numqubits - qubit + 1)
    end
    return circuit
end

"""
    random_circuit_nn(numqubits::Integer, depth::Integer)

Return a random quantum circuit consisting of random unitary gates followed by CNOT or CZ between nearest neighbours. The number of layers is given by depth.

For instance a 2 layer circuit would look like:
-□-┬-□-┬-
-□-☒-□-☒-
-□-┬-□-┬-
-□-☒-□-☒-
with □ random unitaries and -☒ 2 qubits gates like CNOT.
"""
function random_circuit_nn(numqubits::Integer, depth::Integer)::Circuit
    circuit = Circuit(numqubits)
    for layer in 1:depth
        # add 1 qubits random unitary
        for qubit in 1:numqubits
            push!(circuit, GateU(rand(3)...), qubit)
        end
        # add 2 qubit gates, alternate control and target qubits according to parity of layer
        for qubit in (layer%2+1):2:numqubits-1
            push!(circuit, GateCX(), qubit, qubit + 1)
        end
    end
    return circuit
end
