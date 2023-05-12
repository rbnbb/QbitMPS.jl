using Test
using Logging
using QbitMPS
using QbitMPS: mps2statevector, simulate_circuit_imps, simulate_circuit_mpo

filelist = String[
                  # "qft.jl",
                  # "random_circuit.jl",
                  "mpo.jl",
                 ]

include("statevector_checker.jl")

debuglogger = ConsoleLogger(stderr, Logging.Debug)

@testset "QbitMPS" begin
    @testset "$filename" for filename in filelist
    with_logger(debuglogger) do
        @debug "Running $filename"
        st = @timed include(filename)
        @debug "done $(st.time) seconds"
    end  #logger
    end
end

