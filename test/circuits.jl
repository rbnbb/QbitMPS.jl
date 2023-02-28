include("qft.jl")
const numqubits = 10

@testset "Quantum Fourier Transform" for number in rand(10:2^numqubits, 10)
    @test fft_statevector(number, numqubits) ≈ mps_statevector(number, numqubits)
end
