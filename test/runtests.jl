using Test
using QbitMPS
using QbitMPS: mps2statevector, simulate_circuit_imps, simulate_circuit_mpo

filelist = String["qft.jl", "random_circuit.jl", "mpo.jl"]

include("statevector_checker.jl")

@testset "QbitMPS" begin
    @testset "$filename" for filename in filelist
        @debug "Running $filename"
        include(filename)
    end
end
