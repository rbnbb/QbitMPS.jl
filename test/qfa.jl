# Test using QFA properties

function add_with_qfa(a::Integer, b::Integer)
    int2binary(x, n_bits) = Tuple(map(bit -> 1 & (x >> bit), n_bits-1:-1:0))
    nbits = max(a÷2+1, b÷2+1)
    A = int2binary(a, nbits)
    B = int2binary(b, nbits)
    # encode in an initial state
    psi0 = []  # lowest valued bit is first
    for j in 1:nbits
        push!(psi0, "0", string(A[j]), string(B[j]))
    end
    push!(psi0, "0")
    qfa = compile_circuit_to_mpo(generate_circuit(QFA(), 3*nbits+1))
    direct_svd_truncation!(qfa.data)  # compress exactly
    # @info "mpo" qfa
    psi = MPS(siteinds(qfa), (psi0))  # mpo encoding is highest bit first
    apply!(qfa, psi)
    orthogonalize!(psi, 1)
    result = ITensors.sample(psi) .- 1  # ITensors counts from 1
    result = reverse(result[begin:3:end])  # select |C_out> bits
    return sum([result[j]*2^(j-1) for j in 1:nbits+1])
end

@testset "Quantum Full Adder" begin
    tosum = rand(10:20, 2, 10)
    @test 1 == add_with_qfa(0, 1)
    @test 1 == add_with_qfa(1, 0)
    for (i, j) in zip(tosum[1,:], tosum[2,:])
        @test add_with_qfa(i, j) == i + j
    end
end
