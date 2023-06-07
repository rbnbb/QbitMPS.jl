import ITensors.siteinds

include("truncate.jl")

"""
    CircuitMPO

An MPO respresentation of a quantum circuit.
"""
mutable struct CircuitMPO <: ITensors.AbstractMPS
    data::Vector{ITensor}
end

function CircuitMPO(numqubits::Integer)
    qubits = siteinds("Qubit", numqubits)
    # qubits = removetags.(qubits, "Qubit")  # shorter output for debugging
    return CircuitMPO(qubits)
end

"""
    CircuitMPO(numqubits::Integer)

Construct a default CircuitMPO consisting of identity tensors.
"""
function CircuitMPO(qubits)
    mpo = CircuitMPO(ITensor[])
    numqubits = length(qubits)
    l = Vector{Index}(undef, numqubits)
    normcoeff = 1.0 / sqrt(2)
    l[1] = Index(1, ITensors.defaultlinktags(1))
    push!(mpo.data, ITensor(normcoeff * [1 0; 0 1], qubits[1]', qubits[1], l[1]))
    for j in 2:numqubits-1
        l[j] = Index(1, ITensors.defaultlinktags(j))
        push!(mpo.data, normcoeff * delta(qubits[j]', qubits[j]) * delta(l[j], l[j-1]))
    end
    push!(
        mpo.data,
        ITensor(
            normcoeff * [1 0; 0 1],
            qubits[numqubits]',
            qubits[numqubits],
            l[numqubits-1],
        ),
    )
    return mpo
end

siteinds(mpo::CircuitMPO) = [ITensors.siteinds(mpo, j; plev=0)[1] for j in 1:length(mpo)]

function default_maxdim()
    return Inf
end

function svd_kwargs()
    return (cutoff=1e-10,)
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
    Q, S, T = svd(op_2q, (ctrl_idx', ctrl_idx); svd_kwargs()...)
    T = S * T
    return Q, T
end

function _gate2itensors(::GateCX, ctrl_idx::Index, tgt_idx::Index)
    # hard code CX, it is used a lot, first is control and second target:
    # CX = [ I-C  C] * [I X]^T
    l_idx = Index(2; tags="Link")
    return ITensor(Float64[1; 0;; 0; 0;;; 0; 0;; 0; 1], ctrl_idx', ctrl_idx, l_idx),
    ITensor(Float64[1; 0;; 0; 1;;; 0; 1;; 1; 0], tgt_idx', tgt_idx, l_idx)
end

function _gate2itensors(::GateCCX, ctrl_idx1::Index, ctrl_idx2::Index, tgt_idx::Index)
    # hard code CCX
    # CCX = [ I  C]* [I 0; 0 C] * [I X-I]^T
    l_idx = Index(2; tags="Link")
    return ITensor(Float64[1; 0;; 0; 1;;; 0; 0;; 0; 1], ctrl_idx1', ctrl_idx1, l_idx),
    ITensor(
        Float64[1; 0;; 0; 1;;; 0; 0;; 0; 0;;; 0; 0;; 0; 0;;; 0; 0;; 0; 1],
        ctrl_idx2',
        ctrl_idx2,
        l_idx',
        l_idx,
    ),
    ITensor(Float64[1; 0;; 0; 1;;; -1; 1;; 1; -1], tgt_idx', tgt_idx, l_idx')
end

function _add_2_qubit_gate!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    control, target = numqubits .- gate.targets .+ 1  # count qubits from the end
    @assert 0 < control <= numqubits && 0 < target <= numqubits && control != target
    qidx(mpo, n) = (siteinds(mpo, n; plev=1)[1])
    Q, T = _gate2itensors(gate.gate, qidx(mpo, control), qidx(mpo, target))
    l1_idx = inds(Q; tags="Link")[1]
    pm1 = sign(target - control)  # this is +/- 1
    l2_idx = commonind(mpo.data[control], mpo.data[control+pm1])
    C1 = combiner(l1_idx, l2_idx; tags=tags(l2_idx))
    mpo.data[control] = setprime(mpo.data[control] * Q * C1, 1; plev=2)
    for n in control+pm1:pm1:target-pm1
        l2_idx = commonind(mpo.data[n], mpo.data[n+pm1])
        C2 = combiner(l1_idx', l2_idx; tags=tags(l2_idx))
        mpo.data[n] = mpo.data[n] * delta(l1_idx, l1_idx') * C1 * C2
        C1 = C2
        l1_idx = l1_idx'
    end
    mpo.data[target] = setprime(
        mpo.data[target] *
        setprime(T, abs(target - control) - 1; tags="Link") * C1,
        1;
        plev=2,
    )
    return nothing
end

function _add_3_qubit_gate!(mpo::CircuitMPO, gate::CircuitGate)
    numqubits = length(mpo)
    control1, control2, target = numqubits .- gate.targets .+ 1  # count qubits from the end
    @assert 0 < control1 <= numqubits && 0 < target <= numqubits && control1 != target
    qidx(mpo, n) = siteinds(mpo, n; plev=1)[1]
    Q1, Q2, T = _gate2itensors(
        gate.gate,
        qidx(mpo, control1),
        qidx(mpo, control2),
        qidx(mpo, target),
    )
    l1_idx = inds(Q1; tags="Link")[1]
    pm1 = sign(target - control1)  # this is +/- 1
    l2_idx = commonind(mpo.data[control1], mpo.data[control1+pm1])
    C1 = combiner(l1_idx, l2_idx; tags=tags(l2_idx))
    mpo.data[control1] = setprime(mpo.data[control1] * Q1 * C1, 1; plev=2)
    for n in control1+pm1:pm1:target-pm1
        l2_idx = commonind(mpo.data[n], mpo.data[n+pm1])
        C2 = combiner(l1_idx', l2_idx; tags=tags(l2_idx))
        if n != control2
            mpo.data[n] = mpo.data[n] * delta(l1_idx, l1_idx') * C1 * C2
            prime!(Q2; tags="Link")
        else
            # mpo.data[n] = mpo.data[n] * Q2 * C1 * C2
            mpo.data[n] *= setprime(Q2, 1; plev=1, tags="Link")
            mpo.data[n] *= C1
            mpo.data[n] *= C2
            setprime!(mpo.data[n], 1; plev=2)
        end
        # @info "after" mpo.data[n]
        C1 = C2
        l1_idx = l1_idx'
    end
    mpo.data[target] = setprime(
        mpo.data[target] * setprime(T, abs(target - control1) - 1; tags="Link") * C1,
        1;
        plev=2,
    )
    return nothing
end

function apply!(mpo::CircuitMPO, gate::CircuitGate)
    if numqubits(gate) == 1
        _add_single_qubit_gate_to_mpo!(mpo, gate)
    elseif numqubits(gate) == 2
        _add_2_qubit_gate!(mpo, gate)
    else
        _add_3_qubit_gate!(mpo, gate)
    end
    return nothing
end

function compile_circuit_to_mpo!(mpo, circuit::Circuit; maxdim=default_maxdim())
    for gate in circuit
        apply!(mpo, gate)
    end
    if  maxdim < Inf
        direct_svd_truncation!(mpo.data, maxdim)
    end
    return nothing
end

function compile_circuit_to_mpo(circuit::Circuit; maxdim=default_maxdim())
    numqubits = MimiqCircuits.numqubits(circuit)
    mpo = CircuitMPO(numqubits)
    compile_circuit_to_mpo!(mpo, circuit; maxdim)
    return mpo
end

function compile_circuit_to_mpo(circuit::Circuit, qubits; maxdim=default_maxdim())
    mpo = CircuitMPO(qubits)
    compile_circuit_to_mpo!(mpo, circuit; maxdim)
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

function simulator_fidelity(circuit::Circuit; maxdim=default_maxdim())
    mpo::CircuitMPO = compile_circuit_to_mpo(circuit; maxdim)
    qubits = siteinds(mpo)
    psi::MPS = productMPS(qubits, "0")  # Inital state |00...0>
    apply!(mpo, psi; maxdim)
    mpo_inv = compile_circuit_to_mpo(inverse(circuit), qubits; maxdim)
    apply!(mpo_inv, psi; maxdim)
    ket0 = productMPS(qubits, "0")  # Inital state |00...0>
    return real(round(inner(psi, ket0); digits=10))
end

function contract_mpos(O1, O2)
    T = O1[1] * O2[1]
    for j in 2:length(O1)
        T *= O1[j]
        T *= O2[j]
    end
    return T
end

function do_it()
    mpo = CircuitMPO(2)
    _add_2_qubit_gate!(mpo, CircuitGate(GateCX(), 1, 2))
    psi01 = productMPS(siteinds(mpo), ["0", "1"])
    apply!(mpo, psi01)
    # return mps2statevector(psi01)
    return psi01
end
