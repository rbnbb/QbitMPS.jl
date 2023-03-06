using Test
using QbitMPS

filelist = String["qft.jl", "random_circuit.jl"]

@testset "QbitMPS" begin
    @testset "$filename" for filename in filelist
        @debug "Running $filename"
        include(filename)
    end
end
