abstract type QuantumAlgorithm end
struct QFT <: QuantumAlgorithm end
struct QFA <: QuantumAlgorithm end
struct RandomCircuitNN <: QuantumAlgorithm
    depth::Integer
end
struct RandomCircuitLR <: QuantumAlgorithm
    depth::Integer
end

generate_circuit(::QFT, numqubits::Integer) = quantum_fourier_circuit(numqubits)
generate_circuit(::QFA, numqubits::Integer) = quantum_full_adder(numqubits ÷ 3)
generate_circuit(algorithm::RandomCircuitNN, numqubits::Integer) =
    random_circuit_nn(numqubits, algorithm.depth)
generate_circuit(algorithm::RandomCircuitLR, numqubits::Integer) =
    random_circuit_lr(numqubits, algorithm.depth)

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
    circuit = Circuit()
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
─□─┬─□─┬─
─□─☒─□─☒─
─□─┬─□─┬─
─□─☒─□─☒─
with □ random unitaries and -☒ 2 qubits gates like CNOT.
"""
function random_circuit_nn(numqubits::Integer, depth::Integer)::Circuit
    circuit = Circuit()
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

"""
    random_circuit_lr(numqubits::Integer, depth::Integer)

Return a random quantum circuit consisting of random unitary gates followed by CNOT or CZ between randomly selected qubits (i.e. long-range). The number of layers is given by depth.

For instance a one layer circuit could look like:
─□─┬☒☒──
─□─┼┴┼☒─
─□─☒─┴┼─
─□────┴─
with □ random unitaries and -☒ 2 qubits gates like CNOT.
"""
function random_circuit_lr(numqubits::Integer, depth::Integer)::Circuit
    circuit = Circuit()
    for layer in 1:depth
        # add 1 qubits random unitary
        for qubit in 1:numqubits
            push!(circuit, GateU(rand(3)...), qubit)
        end
        # add (potentially) long-range 2 qubit gates
        for target_qubit in 1:numqubits
            control_qubit = rand(1:numqubits-1)
            # make sure control_qubit is not equal to target qubit
            control_qubit = control_qubit < target_qubit ? control_qubit : control_qubit + 1
            push!(circuit, GateCX(), control_qubit, target_qubit)
        end
    end
    return circuit
end

"""
    quantum_full_adder(nadders::Integer)

Return nadders quantum full adder circuits liked.

One quantum full adder has the following circuit:
|Cin> ─────┬─╳─── |S>
|A>  ─┬─┬─│─│─┬─ |A>
|B>  ─┼─╳─┼─┴─╳─ |B>
|0>  ─╳───╳───── |Cout>
where the symbol "─╳" represents a CX/CCX gate.
"""
function quantum_full_adder(nadders::Integer)
    first = 1
    circ = Circuit()
    for j in 1:nadders
        push!(circ, GateCCX(), first + 1, first + 2, first + 3)
        push!(circ, GateCX(), first + 1, first + 2)
        push!(circ, GateCCX(), first, first + 2, first + 3)
        push!(circ, GateCX(), first + 2, first)
        push!(circ, GateCX(), first + 1, first + 2)
        first += 3
    end
    return circ
end
