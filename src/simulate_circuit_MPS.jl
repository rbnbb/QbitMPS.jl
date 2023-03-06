"""
    simulate_circuit(numqubits::Int, circuit::Circuit)

Simulate the application of gates in circuit to an MPS.
Return the final Matrix Product State (MPS).
By, default the initial state is |0...0>, it can be changed with keyword argument.
"""
function simulate_circuit(numqubits::Int, circuit::Circuit; initial_state="0", maxdim=128)::MPS
    qubits = siteinds("Qubit", numqubits)
    psi0::MPS = productMPS(qubits, initial_state)  # Inital state |00...0>
    gate_ops =  _circuit2mpsops(circuit, qubits)
    psi = apply(gate_ops, psi0; maxdim)
    return psi
end

function _circuit2mpsops(circuit::Circuit, qubits)::Vector{ITensor}
    gate2op(gate) = op(qubits, Circuits.matrix(gate), reverse((circuit.nqubits+1) .- gate.targets))
    return map(gate2op, circuit)
end
