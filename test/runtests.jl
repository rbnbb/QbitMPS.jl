using Test
using QbitMPS


filelist = String["circuits.jl"]

@testset "QbitMPS" begin
    @testset "$filename" for filename in filelist
        @debug "Running $filename"
        include(filename)
    end
end
