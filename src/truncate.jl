@enum CompressionType SingleSite TwoSite

function contract_with_dag_self(A::ITensor, neighbour::ITensor)
    return A * dag(prime(A, commonind(A, neighbour)))
end

function is_in_right_canonical_form(M::Vector{ITensor})
    first_is_right_normalized = (M[1]*dag(M[1])).tensor[1] ≈ 1
    if !first_is_right_normalized
        return false
    else
        is_right_normalized(M, j) = contract_with_dag_self(M[j], M[j-1]).tensor ≈ LinearAlgebra.I
        for j in 2:length(M)
            if !is_right_normalized(M, j)
                return false
            end
        end
        return true
    end
end

function is_in_left_canonical_form(M::Vector{ITensor})
    is_left_normalized(M, j) = contract_with_dag_self(M[j], M[j+1]).tensor ≈ I
    for j in 1:length(M)-1
        if !is_left_normalized(M, j)
            return false
        end
    end
    last_is_left_normalized = (M[end]*dag(M[end])).tensor[1] ≈ 1
    return last_is_left_normalized ? true : false
end

function put_in_right_canonical_form!(M::Vector{ITensor})
    if is_in_right_canonical_form(M)
        return nothing
    end
    # sweep right to left and use qr decomposition
    # M[end] = M[end] ./ sqrt((M[end] * dag(M[end])).tensor[1])
    for j in length(M):-1:2
        R, T = qr(M[j], uniqueinds(M[j], M[j-1]); tags=tags(commonind(M[j], M[j-1])))
        M[j] = R
        M[j-1] *= T
    end
    # M[1] = M[1] ./ sqrt((M[1] * dag(M[1])).tensor[1])
end

function direct_svd_truncation!(M::Vector{ITensor}, maxdim=max(maxdim.(M)...))
    if !is_in_right_canonical_form(M)
        put_in_right_canonical_form!(M)
    end
    for j in 1:length(M)-1
        U, S, V = svd(
            M[j],
            uniqueinds(M[j], M[j+1]);
            lefttags=tags(commonind(M[j], M[j+1])),
            cutoff=1e-18,
            maxdim,
        )
        M[j] = U
        M[j+1] *= S * V
    end
end

"""
Move orthogonality center from site j one site to left or right, as indicated by step.
"""
function move_orthogonality_center!(M::Vector{ITensor}, j::Integer, step::Integer)
    @assert step == 1 || step == -1 "This function only moves center by one site."
    A, B = qr(M[j], uniqueinds(M[j], M[j+step]); tags=tags(commonind(M[j], M[j+step])))
    M[j] = A
    return M[j+step] *= B
end

"""
Compute bipartite entropy, provided T is orthogonality center
"""
function von_neumann(T::ITensor, inds)
    _, S, _ = svd(T, inds)
    return sum([x != 0 ? -2 * x^2 * log2(x) : 0 for x in diag(S)])
end

function contractions_for_sweep(
    psi::Vector{ITensor},
    phi::Vector{ITensor},
    range,
)::Vector{ITensor}
    R = Vector{ITensor}(undef, length(range))
    R[last(range)-step(range)] = psi[last(range)] * phi[last(range)]
    for j in reverse(range[begin:end-2])
        R[j] = R[j+step(range)] * phi[j+step(range)] * psi[j+step(range)]
    end
    return R
end

"""
Variationnally compress state psi to maxdim, starting from guess state phi.
The site indices of phi need to be the same as those of psi.
"""
function variationally_compress(
    psi::Vector{ITensor},
    maxdim::Integer,
    nsweeps::Integer=1,
    algo::CompressionType=TwoSite;
    fidelities=nothing,
    entropies=nothing,
)
    n_points_around(x, n) = x-div(n, 2)-n%2+1:x+div(n, 2)
    """ Sweep once through MPS using single site variational compression algorithm.
    B is a compressed version of A. Optionally return information such as fidelity and entanglement entropy. """
    function sweep_once_1s!(
        range,
        A::Vector{ITensor},
        B::Vector{ITensor};
        compute_fidelity=false,
        compute_entropy=false,
    )
        beg = first(range)  # for readability
        overlap_A_B, entropy = missing, missing
        B .= dag.(B)  # complex conjugate guess state, per algorithm
        R = contractions_for_sweep(A, B, range)
        if compute_fidelity
            # compute overlap BEFORE sweep using R, it's cheap
            overlap_A_B = real((R[beg]*A[beg]*B[beg]).tensor[1])
        end
        if compute_entropy  # average entropy over multiple central bonds
            avg_range = n_points_around(length(range) ÷ 2, 3)
            entropy = 0.0
        end
        # Note on the use of dag(): by conjugating again each result
        # we avoid unnecessary assignments and save memory
        B[beg] = dag(A[beg] * R[beg])
        move_orthogonality_center!(B, beg, step(range))
        L = A[beg] * B[beg]
        for j in range[begin+1:end-1]
            B[j] = dag(L * A[j] * R[j])
            # compute von Neumann entropy of middle bipartition
            compute_entropy && j in avg_range &&
                (entropy += von_neumann(B[j], commonind(B[j], B[j+1])))
            move_orthogonality_center!(B, j, step(range))
            L *= B[j]
            L *= A[j]
        end
        B[last(range)] = dag(L * A[last(range)])
        B .= dag.(B)  # put state in normal mode
        return overlap_A_B, entropy / length(avg_range)
    end
    function sweep_once_2s!(
        range,
        A::Vector{ITensor},
        B::Vector{ITensor};
        compute_fidelity=false,
        compute_entropy=false,
    )
        beg = first(range)  # for readability
        one = step(range)  # for readability
        overlap_A_B, entropy = missing, missing
        B .= dag.(B)  # complex conjugate guess state, per algorithm
        R = contractions_for_sweep(A, B, range)
        if compute_fidelity
            # compute overlap BEFORE sweep using R, it's cheap
            overlap_A_B = real((R[beg]*A[beg]*B[beg]).tensor[1])
        end
        if compute_entropy  # average entropy over multiple central bonds
            avg_range = n_points_around(length(range) ÷ 2, 3)
            avg_range = avg_range .+ 1
            entropy = 0.0
        end
        T = A[beg] * A[beg+one] * R[beg+one]
        B[beg], _ = svd(T,
            uniqueinds(T, B[beg+one]);
            lefttags=tags(commonind(A[beg], A[beg+one])), maxdim)
        B[beg] = dag(B[beg])
        L = A[beg] * B[beg]
        for j in range[begin+1:end-2]
            T = (L * A[j]) * (A[j+one] * R[j+one])
            B[j], _ =
                svd(T,
                    uniqueinds(T, B[j+one]);
                    lefttags=tags(commonind(B[j], B[j+one])), maxdim)
            B[j] = dag(B[j])
            # compute von Neumann entropy of middle bipartition
            compute_entropy && j in avg_range &&
                (entropy += von_neumann(B[j], commonind(B[j], B[j-one])))
            L *= B[j]
            L *= A[j]
        end
        T = (L * A[last(range)-one]) * A[last(range)]
        B[last(range)-one], S, V = svd(T,
            uniqueinds(T, B[last(range)]);
            lefttags=tags(commonind(A[last(range)-one], A[last(range)])), maxdim)
        B[last(range)-one] = dag(B[last(range)-one])
        # L *= A[last(range)-one]
        # L *= B[last(range)-one]
        # B[last(range)] = dag(L * A[last(range)])
        B[last(range)] = dag(S * V)
        B .= dag.(B)  # put state in normal mode
        return overlap_A_B, entropy / length(avg_range)
    end
    # @info "$(string(algo)) var compress from $(max(ITensors.maxdim.(psi)...)) to $maxdim \n"
    sweep_once! = algo == SingleSite ? sweep_once_1s! : sweep_once_2s!
    phi = deepcopy(psi)
    direct_svd_truncation!(phi, maxdim)
    put_in_right_canonical_form!(phi)
    left2right = 1:length(psi)
    right2left = length(psi):-1:1
    F = -42.0
    S = -42.0
    for j in 1:nsweeps
        # @debug sv(psi)' * sv(phi)
        F, S = sweep_once!(
            isodd(j) ? left2right : right2left,
            psi,
            phi;
            compute_entropy=!isnothing(entropy),
            compute_fidelity=!isnothing(fidelities),
        )
        @debug F
        @debug S
        isnothing(fidelities) || (fidelities[j] = F)
        isnothing(entropies) || (entropies[j] = S)
    end
    return phi
end

function fidelity_entropy_per_sweep(
    numqubits::Integer,
    nsweeps::Integer,
    state="random",
    compress_to=nothing,
    algo=SingleSite,
)
    fidelities = Vector{Float64}(undef, nsweeps)
    entropies = Vector{Float64}(undef, nsweeps)
    if state == "random"
        χ = 2^(numqubits ÷ 2)
        psi = randomMPS(ComplexF64, siteinds("Qubit", numqubits); linkdims=χ)
        direct_svd_truncation!(psi.data, maxlinkdim(psi))
    elseif isinteger(state)  # "shallow_circuit with depth state"
        psi = simulate_circuit(
            generate_circuit(RandomCircuitNN(state), numqubits),
            QbitMPS.sMPS,
        )
    end
    isnothing(compress_to) && (compress_to = maxlinkdim(psi) ÷ 3 + 1)
    variationally_compress(psi.data, compress_to, nsweeps, algo; fidelities, entropies)
    return (fidelities=fidelities, entropies=entropies)
end

#=
include("utils.jl")
sv = mps2statevector
@testset "Orthogonality tools" begin
    function rnd_mps(N, ortho)
        psi = randomMPS(ComplexF64, siteinds("Qubit", N); linkdims=500)
        orthogonalize!(psi, ortho)
        return psi.data
    end
    @test is_in_left_canonical_form(rnd_mps(5, 5))
    @test is_in_right_canonical_form(rnd_mps(5, 1))
end
=#
