using Test
using QbitMPS

filelist = String["qft.jl", "random_circuit.jl", "mpo.jl"]

include("statevector_checker.jl")

@testset "QbitMPS" begin
    @testset "$filename" for filename in filelist
        @debug "Running $filename"
        include(filename)
    end
end
