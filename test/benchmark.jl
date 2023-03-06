using RecipesBase
using JLD2

push!(LOAD_PATH, pwd())
using QbitMPS

include("qft.jl")

@enum mps_method imps
@enum check_method FFT statevector

const file_storage_id = 1

"""
    OneRun

Wrapper for data about one calculation and comparison to exact method.

# Fields

  - `numqubits::Int`: number of qubits
  - `max_bond_dimension::Int`: parameter for MPS simulation
  - `sigma::Float64`: standard deviation between MPS and exact simulation statevectors
  - `time::Float64`: MPS computation duration in seconds
  - `bytes::Int`: memory used for MPS simulation
  - `check_time::Float64`: statevector computation duration in seconds
  - `check_bytes::Int`: memory used for statevector simulation
"""
struct OneBenchmark
    numqubits::Int
    max_bond_dimension::Int
    sigma::Float64
    fidelity::Float64
    depth::Union{Int,Missing}  # for circuits of variable size
    time::Float64
    bytes::Int
    check_time::Float64
    check_bytes::Int
end

"""
    Results

Wrapper for data to facilitate plotting and storage.

# Fields

  - `algorithm::quantum_algorithm`:
  - `method::mps_method`:
  - `checker::check_method`:
  - `runs::Vector{OneRun}`:
"""
struct BenchmarkResults
    algorithm::Type{<:QuantumAlgorithm}
    method::mps_method
    checker::check_method
    runs::Vector{OneBenchmark}
end

function run_benchmark(
    algorithm::Type{<:QuantumAlgorithm},
    numqubits,
    max_bond_dimension,
    depth = missing,
)::BenchmarkResults
    runs = _run_one_benchmark.(algorithm, numqubits, max_bond_dimension, depth)
    r = BenchmarkResults(algorithm, imps, FFT, runs)
    result_name = "$(string(algorithm))_n=$(numqubits)_chi=$(max_bond_dimension)"
    fname = "./test/output/data_$(file_storage_id).jld2"
    jldopen(fname, "a+") do file
        file[result_name] = r
    end
    return r
end

function open_datafile()
    return jldopen("./test/output/data_$(file_storage_id).jld2", "r")
end

function _run_one_benchmark(
    ::Type{QFT},
    numqubits::Integer,
    max_bond_dimension::Integer,
    depth = missing,
)
numstate = rand(0:2^numqubits)
mps_stats = @timed mps_sv = mps_statevector(numstate, numqubits; max_bond_dimension)
    fft_stats = @timed fft_sv = fft_statevector(numstate, numqubits)
    sigma = sqrt(abs(sum((mps_sv - fft_sv) .^ 2) / (2^numqubits)))
    fidelity = abs((fft_sv' * mps_sv))^2
    return OneBenchmark(
        numqubits,
        max_bond_dimension,
        sigma,
        fidelity,
        missing,
        mps_stats.time,
        mps_stats.bytes,
        fft_stats.time,
        fft_stats.bytes,
    )
end

@recipe f(::Type{BenchmarkResults}, r::BenchmarkResults) = r.runs

@recipe function f(results::BenchmarkResults, xaxis, yaxis)
    xlabel --> string(xaxis)
    ylabel --> string(yaxis)
    label --> string(results.algorithm) * " " * string(yaxis)
    markershape --> :circle
    xs = getproperty.(results.runs, (xaxis))
    ys = getproperty.(results.runs, (yaxis))
    xs, ys
end
