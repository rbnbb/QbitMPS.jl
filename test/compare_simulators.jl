using QbitMPS: von_neumann, simulator_fidelity
using DataFrames
using CSV
using Random

"""Compute bipartite entropy by averaging over a few central bonds"""
function entropy(psi::MPS)
    n_points_around(x, n) = x-div(n, 2)-n%2+1:x+div(n, 2)
    n_avg = 3
    S = 0.0
    for j in n_points_around(length(psi)÷2, n_avg)
        orthogonalize!(psi, j)
        S += von_neumann(psi[j], commonind(psi[j], psi[j+1]))
    end
    return round(S/n_avg; digits=10)
end

function entanglement_x_depth(N, depths, chis)
    # N = 15
    entropies = Float64[]
    df = DataFrame(D=Int[], χ=Int[], S=Float64[])
    for chi in chis
        for D in depths
            Random.seed!(1)
            c = generate_circuit(RandomCircuitNN(D), N)
            mpo = compile_circuit_to_mpo(c; maxdim=chi)
            psi = MPS(siteinds(mpo), "0")
            apply!(mpo, psi)
            push!(df, (D, chi, entropy(psi)))
        end
    end
    CSV.write("test/output/S_x_D$(depths)_4_Chi$(chis).csv", df)
    return df
end

function fidelity_x_N(Ns, chis, D)
    df = DataFrame(N=Int[], χ=Int[], F=Float64[], t=Float64[])
    # D = 10
    for chi in chis
        for N in Ns
            Random.seed!(1)
            c = generate_circuit(RandomCircuitNN(D), N)
            stats = @timed fidelity = simulator_fidelity(c; maxdim=chi)
            push!(df, (N, chi, fidelity, stats.time))
        end
    end
    # end
    CSV.write("test/output/F_x_N$(Ns)_4_Chi$(chis).csv", df)
    return df
end

# function time_x_N_sv(Ns, D)
#     df = DataFrame(N=Int[], tsv=Float64[]
#     for n in 1:max(Ns..., 20)
#             Random.seed!(1)
#             c = generate_circuit(RandomCircuitNN(D), N)
#             stats = @timed 
#     end

function S_x_D(df)
    chis = sort(collect(Set(df.χ)))
    Ds = df[df.χ .== chis[1], :D]
    S(chi) = df[df.χ .== chi, :S]
    kwargs = (
             xlab = "# layers",
             ylab = "S",
             lab = reshape(["χ = $(chi)" for chi in chis], (1, length(chis))),
             )
    plot(Ds, [S(chi) for chi in chis]; kwargs...)
end

function F_x_N(df)
    chis = sort(collect(Set(df.χ)))
    Ns = df[df.χ .== chis[1], :N]
    F(chi) = df[df.χ .== chi, :F]
    kwargs = (
             xlab = "# qubits",
             ylab = "fidelity",
             lab = reshape(["χ = $(chi)" for chi in chis], (1, length(chis))),
             markershape = [:utriangle :dtriangle :rtriangle],
             # linestyle = :dot,
             ms = 3,
             ma = 0.9,
             ylim =  (-0.02, 1.04)
             )
    plot(Ns, [F(chi) for chi in chis]; kwargs...)
end

plot_csv_df(input, plotter) = plotter(DataFrame(CSV.File(input)))
