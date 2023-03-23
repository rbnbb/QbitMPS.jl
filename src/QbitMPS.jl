"""
QbitMPS is a project that tries to explore simulating quantum computing using
classical Matrix Product State (MPS) methods.
"""
module QbitMPS

using ITensors
using Circuits

include("utils.jl")
include("circuits.jl")
include("imps.jl")
include("mpo.jl")

export QuantumAlgorithm, QFT, RandomCircuitNN

export generate_circuit
export simulate_circuit

end # module QbitMPS
