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
    push!(mpo.data, ITensor([1 0; 0 1], qubits[1]', qubits[1], l[1]))
    for j in 2:numqubits-1
        l[j] = Index(1, ITensors.defaultlinktags(j))
        push!(mpo.data, delta(qubits[j]', qubits[j]) * delta(l[j], l[j-1]))
    end
    push!(
        mpo.data,
        ITensor([1 0; 0 1], qubits[numqubits]', qubits[numqubits], l[numqubits-1]),
    )
    return mpo
end

siteinds(mpo::CircuitMPO) = [ITensors.siteinds(mpo, j; plev = 0)[1] for j in 1:length(mpo)]

function _add_single_qubit_gate_to_mpo!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    target = numqubits - gate.targets[1] + 1  # count qubits from the end
    site = siteinds(mpo, target; plev = 1)[1]
    op_1q = ITensor(Circuits.matrix(gate.gate), site', site)  # BE CAREFUL where  you put the prime
    mpo.data[target] = setprime(mpo.data[target] * op_1q, 1; plev = 2)
    return nothing
end

function _add_2_qubit_gate!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    control, target = numqubits .- gate.targets .+ 1  # count qubits from the end
    @assert 0 < control <= numqubits && 0 < target <= numqubits && control != target
    ctrl_idx = siteinds(mpo, control; plev = 1)[1]
    tgt_idx = siteinds(mpo, target; plev = 1)[1]
    op_2q = ITensor(Circuits.matrix(gate.gate), tgt_idx', ctrl_idx', tgt_idx, ctrl_idx)
    # Q, T = qr(op_2q, (ctrl_idx', ctrl_idx))
    Q, S, T = svd(op_2q, (ctrl_idx', ctrl_idx); cutoff = 1e-18)
    T = S * T
    l1_idx = inds(Q; tags = "Link")[1]
    pm1 = sign(target - control)  # this is +/- 1
    l2_idx = commonind(mpo.data[control], mpo.data[control+pm1])
    C1 = combiner(l1_idx, l2_idx; tags = tags(l2_idx))
    mpo.data[control] = setprime(mpo.data[control] * Q * C1, 1; plev = 2)
    for n in control+pm1:pm1:target-pm1
        l2_idx = commonind(mpo.data[n], mpo.data[n+pm1])
        C2 = combiner(l1_idx', l2_idx; tags = tags(l2_idx))
        mpo.data[n] = mpo.data[n] * delta(l1_idx, l1_idx') * C1 * C2
        C1 = C2
        l1_idx = l1_idx'
    end
    mpo.data[target] = setprime(
        mpo.data[target] * setprime(T, abs(target - control) - 1; tags = "Link") * C1,
        1;
        plev = 2,
    )
end

function do_it()
    mpo = CircuitMPO(2)
    _add_2_qubit_gate!(mpo, CircuitGate(GateCX(), 1, 2))
    psi01 = productMPS(siteinds(mpo), ["0", "1"])
    apply!(mpo, psi01)
    # return mps2statevector(psi01)
    return psi01
end

function apply!(mpo::CircuitMPO, gate::CircuitGate)
    if numqubits(gate) == 1
        _add_single_qubit_gate_to_mpo!(mpo, gate)
    elseif numqubits(gate) == 2
        _add_2_qubit_gate!(mpo, gate)
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
    psi.data[1] = noprime(psi.data[1] * mpo.data[1] * C1)
    N = length(psi)
    for j in 2:N-1
        C2 = C1
        C1 = combiner(l_psi[j], l_mpo[j]; tags = tags(l_psi[j]))
        psi.data[j] = noprime(psi.data[j] * mpo.data[j] * C2 * C1)
    end
    psi.data[N] = noprime(psi.data[N] * mpo.data[N] * C1)
    return nothing
end

function simulate_circuit_mpo(circuit::Circuit, maxdim = 128)
    mpo::CircuitMPO = compile_circuit_to_mpo(circuit)
    qubits = siteinds(mpo)
    psi0::MPS = productMPS(qubits, "0")  # Inital state |00...0>
    apply!(mpo, psi0; maxdim)
    return psi0
end

function random_circuit_only_1q_unitaries(numqubits::Integer)
    circuit = Circuit(numqubits)
    for qubit in 1:numqubits
        push!(circuit, GateU(rand(3)...), qubit)
    end
    return circuit
end

@testset "Test 1 qubit gates for MPO simulator" begin
    function test_1qubit_unitaries(numqubits::Integer)
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
        circ = Circuit(2)
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
    MAXQUBITS = 7
    CIRCUIT_DEPTH = 2
    for n in 2:MAXQUBITS
        @test test_2q_nn(n, CIRCUIT_DEPTH)
    end
end

@testset "Test 2 qubits long-range gates for MPO simulator" begin
    function test_bell_preparation_3q()
        circ = Circuit(3)
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
        circ = random_circuit(numqubits, depth)
        mps_sv = mps2statevector(simulate_circuit_mpo(circ))
        sv = check_statevector(circ, max(4, numqubits))[1:2^numqubits]
        return mps_sv ≈ sv
    end
    @test test_bell_preparation_3q()
    MAXQUBITS = 7
    CIRCUIT_DEPTH = 2
    for n in 2:MAXQUBITS
        @test test_2q_long_range(n, CIRCUIT_DEPTH)
    end
end
