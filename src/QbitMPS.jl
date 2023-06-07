"""
QbitMPS is a project that tries to explore simulating quantum computing using
classical Matrix Product State (MPS) methods.
"""
module QbitMPS

using ITensors
using MimiqCircuits

include("utils.jl")
include("circuits.jl")
include("imps.jl")
include("mpo.jl")
include("simulator.jl")

export QuantumAlgorithm, QFT, QFA, RandomCircuitNN, RandomCircuitLR

export generate_circuit
export Simulator, sMPO, sMPS
export simulate_circuit
export compile_circuit_to_mpo
export apply!

end # module QbitMPS
