import ITensors.siteinds

include("truncate.jl")

"""
    CircuitMPO

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
    normcoeff = 1.0 / sqrt(2)
    l[1] = Index(1, ITensors.defaultlinktags(1))
    push!(mpo.data, ITensor( normcoeff*[1 0; 0 1], qubits[1]', qubits[1], l[1]))
    for j in 2:numqubits-1
        l[j] = Index(1, ITensors.defaultlinktags(j))
        push!(mpo.data,  normcoeff*delta(qubits[j]', qubits[j]) * delta(l[j], l[j-1]))
    end
    push!(
        mpo.data,
        ITensor( normcoeff*[1 0; 0 1], qubits[numqubits]', qubits[numqubits], l[numqubits-1]),
    )
    return mpo
end

siteinds(mpo::CircuitMPO) = [ITensors.siteinds(mpo, j; plev=0)[1] for j in 1:length(mpo)]

function default_maxdim()
    return 64
end

function _add_single_qubit_gate_to_mpo!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    target = numqubits - gate.targets[1] + 1  # count qubits from the end
    site = siteinds(mpo, target; plev=1)[1]
    op_1q = ITensor(MimiqCircuits.matrix(gate.gate), site', site)  # BE CAREFUL where  you put the prime
    mpo.data[target] = setprime(mpo.data[target] * op_1q, 1; plev=2)
    return nothing
end

function _gate2itensors(gate::AbstractGate{2}, ctrl_idx::Index, tgt_idx::Index)
    op_2q = ITensor(MimiqCircuits.matrix(gate), tgt_idx', ctrl_idx', tgt_idx, ctrl_idx)
    # Q, T = qr(op_2q, (ctrl_idx', ctrl_idx))
    Q, S, T = svd(op_2q, (ctrl_idx', ctrl_idx); cutoff=1e-10)
    T = S * T
    return Q, T
end

function _gate2itensors(gate::GateCX, ctrl_idx::Index, tgt_idx::Index)
    # hard code CX, it is used a lot, first is control and second target:
    # CX = [ I-C  C] * [I X]^T
    l_idx = Index(2; tags="Link")
    return ITensor(Float64[1;0;;0;0;;;0;0;;0;1], ctrl_idx', ctrl_idx, l_idx),ITensor(Float64[1;0;;0;1;;;0;1;;1;0], tgt_idx', tgt_idx, l_idx)
end

function _gate2itensors(gate::GateCCX, ctrl_idx1::Index, ctrl_idx2::Index, tgt_idx::Index)
    # hard code CCX
    # CCX = [ I-C  C]* [I 0; 0 C] * [I X]^T
    l_idx = Index(2; tags="Link")
    return ITensor(Float64[1;0;;0;0;;;0;0;;0;1], ctrl_idx1', ctrl_idx1, l_idx), ITensor(Float64[1;0;;0;1;;;0;0;;0;0;;;0;0;;0;0;;;0;0;;0;1], ctrl_idx2', ctrl_idx2, l_idx), ITensor(Float64[1;0;;0;1;;;0;1;;1;0], tgt_idx', tgt_idx, l_idx)
end

function _add_2_qubit_gate!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    normcoeff = 1.0 #/ sqrt(2)
    control, target = numqubits .- gate.targets .+ 1  # count qubits from the end
    @assert 0 < control <= numqubits && 0 < target <= numqubits && control != target
    ctrl_idx = siteinds(mpo, control; plev=1)[1]
    tgt_idx = siteinds(mpo, target; plev=1)[1]
    Q, T = _gate2itensors(gate.gate, ctrl_idx, tgt_idx)
    l1_idx = inds(Q; tags="Link")[1]
    pm1 = sign(target - control)  # this is +/- 1
    l2_idx = commonind(mpo.data[control], mpo.data[control+pm1])
    C1 = combiner(l1_idx, l2_idx; tags=tags(l2_idx))
    mpo.data[control] = setprime(normcoeff * mpo.data[control] * Q * C1, 1; plev=2)
    for n in control+pm1:pm1:target-pm1
        l2_idx = commonind(mpo.data[n], mpo.data[n+pm1])
        C2 = combiner(l1_idx', l2_idx; tags=tags(l2_idx))
        mpo.data[n] = normcoeff * mpo.data[n] * delta(l1_idx, l1_idx') * C1 * C2
        C1 = C2
        l1_idx = l1_idx'
    end
    mpo.data[target] = setprime(
        normcoeff * mpo.data[target] * setprime(T, abs(target - control) - 1; tags="Link") * C1,
        1;
        plev=2,
    )
    return nothing
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

# function truncate_mpo_svd!(mpo::CircuitMPO)
#     nothing
# end
function compile_circuit_to_mpo!(mpo, circuit::Circuit; maxdim=default_maxdim())
    for gate in circuit
        apply!(mpo, gate)
        # if max(linkdims(mpo)...) > maxdim
        #     truncate_mpo_svd!(mpo)
        # end
    end
    return nothing
end

function compile_circuit_to_mpo(circuit::Circuit; maxdim=default_maxdim())
    numqubits = MimiqCircuits.numqubits(circuit)
    mpo = CircuitMPO(numqubits)
    for gate in circuit
        apply!(mpo, gate)
        # if max(linkdims(mpo)...) > maxdim
        #     truncate_mpo_svd!(mpo)
        # end
    end
    return mpo
end

function apply!(mpo::CircuitMPO, psi::MPS; kwargs...)
    l_psi = linkinds(psi)
    l_mpo = linkinds(mpo)
    C1 = combiner(l_psi[1], l_mpo[1]; tags=tags(l_psi[1]))
    normcoeff = sqrt(2)
    psi.data[1] = noprime(psi.data[1] * mpo.data[1] * normcoeff * C1)
    N = length(psi)
    for j in 2:N-1
        C2 = C1
        C1 = combiner(l_psi[j], l_mpo[j]; tags=tags(l_psi[j]))
        psi.data[j] = noprime(psi.data[j] * mpo.data[j] * normcoeff * C2 * C1)
    end
    psi.data[N] = noprime(psi.data[N] * mpo.data[N] * normcoeff * C1)
    return nothing
end

function simulate_circuit_mpo(circuit::Circuit; maxdim=default_maxdim())
    mpo::CircuitMPO = compile_circuit_to_mpo(circuit; maxdim)
    qubits = siteinds(mpo)
    psi0::MPS = productMPS(qubits, "0")  # Inital state |00...0>
    apply!(mpo, psi0; maxdim)
    return psi0
end

function contract_mpos(O1, O2)
    T = O1[1] * O2[1]
    for j in 2:length(O1)
        T *= O1[j]
        T *= O2[j]
    end
    return T
end
