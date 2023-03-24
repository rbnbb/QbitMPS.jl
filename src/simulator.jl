@enum Simulator sMPS sMPO

"""
    simulate_circuit(circuit::Circuit, method::Simulator; kwargs...)

Return an ITensors::MPS with the quantum state after the execution of the circuit.
Number of qubits is take from circuit.
"""
function simulate_circuit(circuit::Circuit, method::Simulator; kwargs...)::MPS
    if method == sMPS
        return simulate_circuit_imps(circuit; kwargs...)
    elseif method == sMPO
        return simulate_circuit_mpo(circuit; kwargs...)
    end
    throw(
        ArgumentError(
            "Invalid method provided. Please check available methods with instances(Simulator)",
        ),
    )
end
