using ITensors
using Circuits

import ITensors.siteinds
include("circuits.jl")
include("utils.jl")
include("../test/statevector_checker.jl")

"""
    CircMPO

An MPO respresentation of a quantum circuit.
"""
mutable struct CircuitMPO <: ITensors.AbstractMPS
    data::Vector{ITensor}
end

"""
    CircuitMPO(numqubits::Integer)

Construct a default CircuitMPO consisting of identity tensors.
"""
function CircuitMPO(numqubits::Integer)
    mpo = CircuitMPO(ITensor[])
    qubits = siteinds("Qubit", numqubits)
    l = Vector{Index}(undef, numqubits)
    l[1] = Index(1, ITensors.defaultlinktags(1))
    push!(mpo.data, delta(qubits[1]', qubits[1], l[1]))
    for j in 2:numqubits-1
        l[j] = Index(1, ITensors.defaultlinktags(j))
        push!(mpo.data, delta(qubits[j]', qubits[j], l[j], l[j-1]))
    end
    push!(mpo.data, delta(qubits[numqubits]', qubits[numqubits], l[numqubits-1]))
    return mpo
end

siteinds(mpo::CircuitMPO) = [ITensors.siteinds(mpo, j; plev = 0)[1] for j in 1:length(mpo)]

function _add_single_qubit_gate_to_mpo!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    target = numqubits - gate.targets[1] + 1  # count qubits from the end
    site = siteinds(mpo, target; plev = 1)[1]
    op = ITensor(Circuits.matrix(gate.gate), prime(site), site)  # BE CAREFUL where  you put the prime
    mpo.data[target] = setprime(mpo.data[target] * op, 1; plev = 2)
    return nothing
end

function apply!(mpo::CircuitMPO, gate::CircuitGate)
    if numqubits(gate) == 1
        _add_single_qubit_gate_to_mpo!(mpo, gate)
    elseif numqubits(gate == 2)
        return nothing
    end
    return nothing
end

function compile_circuit_to_mpo(circuit::Circuit)
    numqubits = circuit.nqubits
    mpo = CircuitMPO(numqubits)
    for gate in circuit
        apply!(mpo, gate)
    end
    return mpo
end

function apply!(mpo::CircuitMPO, psi::MPS; kwargs...)
    l_psi = linkinds(psi)
    l_mpo = linkinds(mpo)
    C1 = combiner(l_psi[1], l_mpo[1]; tags = tags(l_psi[1]))
    psi[1] = noprime(psi[1] * mpo[1] * C1)
    N = length(psi)
    for j in 2:N-1
        C2 = C1
        C1 = combiner(l_psi[j], l_mpo[j]; tags = tags(l_psi[j]))
        psi[j] = noprime(psi[j] * mpo[j] * C2 * C1)
    end
    psi[N] = noprime(psi[N] * mpo[N] * C1)
    return nothing
end

function simulate_circuit_mpo(circuit::Circuit, maxdim = 128)
    mpo::CircuitMPO = compile_circuit_to_mpo(circuit)
    qubits = siteinds(mpo)
    psi0::MPS = productMPS(qubits, "0")  # Inital state |00...0>
    apply!(mpo, psi0; maxdim)
    return psi0
end

function random_circuit_only_unitaries(numqubits::Integer)
    circuit = Circuit(numqubits)
    for qubit in 1:numqubits
        push!(circuit, GateU(rand(3)...), qubit)
    end
    return circuit
end

function test_1qubit_unitaries(numqubits::Integer)
    @assert numqubits > 1
    circ = random_circuit_only_unitaries(numqubits)
    if numqubits == 2
        sv = check_statevector(circ, 4)[1:2^numqubits]
    elseif numqubits == 3
        sv = check_statevector(circ, 4)[1:2^numqubits]
    else
        sv = check_statevector(circ, numqubits)
    end
    mps_sv = mps2statevector(simulate_circuit_mpo(circ))

    stringify(ind) = [string(string.(j...)) for j in ind]

    # return [stringify(list_fock_states(numqubits)) round.(sv; digits = 3) round.( mps_sv; digits = 3,)]
    return mps_sv ≈ sv
end

@testset "Random unitaries circuit test for MPO simulator" begin
    MAXQUBITS = 7
    for n in 3:MAXQUBITS
        @test test_1qubit_unitaries.(n)
    end
end
