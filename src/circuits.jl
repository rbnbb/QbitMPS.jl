@doc raw"""
    quantum_fourier_circuit(numqubits::Int)

Return the Quantum Fourier Transform circuit for numqubits.

The formula for QFT used is:
```math
QFT_N | x \rangle = \frac{1}{\sqrt{N}} \sum_{y=0}^{N-1} e^{\frac{2 \pi i xy }{ 2^n}} | y \rangle
```
with ``|x\rangle`` a number state and ``N=2^{\text{numqubits}}``.
"""
function quantum_fourier_circuit(numqubits::Int)::Circuit
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
