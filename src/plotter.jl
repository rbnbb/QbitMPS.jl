# Hold functions for plotting here for reproducibility
# Use the naming scheme
# [multi_]xaxis_yaxis[4multi]
using Plots

function fidelity_sweep_1s_2s(numqubits, nsweeps)
    Random.seed!(1)
    f_1s, _ = fidelity_entropy_per_sweep(numqubits, nsweeps, "random", 16, SingleSite)
    Random.seed!(1)
    f_2s, _ = fidelity_entropy_per_sweep(numqubits, nsweeps, "random", 16, TwoSite)
    kwargs = (
        # ylims=(0, 1),
        title="Variational compression from $(2^(numqubits÷2)) to 16\n N=$numqubits, $nsweeps sweeps",
        xlabel="# sweep",
        ylabel="fidelity",
    )
    return plot([f_1s, f_2s]; label=["1 site" "2 site"], kwargs...)
end

function fgain4diffstates_sweep(numqubits, nsweeps, nlines)
    # run a lot of random states with :
    Random.seed!(1)
    fidelities = []
    for i in 1:nlines
        f, _ = fidelity_per_sweep(numqubits, nsweeps)
        push!(fidelities, f[2:end] .- f[1:end-1])
    end
    return plot(
        fidelities;
        c=:viridis,
        line_z=(1:nlines)',
        yscale=:log10,
        title="Single site variational compression\n $nlines randomMPS; 50 sweeps; N=12",
        legend=false,
        xlabel="# sweep",
        ylabel="fidelity gain",
    )
end

function fgain4D_sweep(numqubits::Integer, nsweeps::Integer, depth=[5, 8, 16, 32])
    ys = Matrix{Float64}(undef, nsweeps, length(depth))
    label = Matrix{String}(undef, 1, length(depth))
    for (j, D) in enumerate(depth)
        Random.seed!(3)
        n_avg = 5
        fidelities = sum([fidelity_entropy_per_sweep(numqubits, nsweeps, D)[:fidelities] for j in 1:n_avg])./n_avg
        # @show fidelities = round.(fidelities; digits = 5)
        # if all([fidelities[j] > fidelities[j-1] for j in 2:length(fidelities)])
        #     print("YES")
        #     ys[:,j] = vcat([42.0], fidelities[2:end] .- fidelities[1:end-1])
        # else
        #     # ys[:,j] = fidelities
        # end
        ys[:,j] = fidelities
        label[j] = "D=$D"
    end
    f_kwargs = (
        xlabel="# sweep",
        ylabel="fidelity",# gain",
        markershape=:circle,
        markersize=1,
        # yscale=:log10,
        # yticks=[10.0^j for j in -11:2:3],
        formatter=:plain,
        title="Variational compression",
    )
    return plot(ys; label, f_kwargs...)
end

function F_maxdim(numqubits::Integer)
    nsweeps = 5
    nstates2avg = 10
    maxdims = 2:4:2^(numqubits÷2)
    fidelities = Vector{Float64}(undef, length(maxdims))
    entropies = Vector{Float64}(undef, length(maxdims))
    for (j, χ) in enumerate(maxdims)
        Random.seed!(1)
        fidelities[j], entropies[j] = 0, 0
        for _ in 1:nstates2avg
            results = fidelity_per_sweep(numqubits, nsweeps, "random", χ)
            fidelities[j] += results[:fidelities][end]
            entropies[j] += median(results[:entropies]) / log2(numqubits)
        end
        fidelities[j] /= nstates2avg
        entropies[j] /= nstates2avg
    end
    kwargs = (dpi=400, xflip=true, formatter=:plain, ms=1, legend=false)
    f_kwargs =
        (xlabel="χ", ylabel="fidelity", label="fidelity",
            title="Variational compression, 1 site\n $nstates2avg randomMPS/point; $nsweeps sweeps; N=$(numqubits)",
            markershape=:circle, color=:red)
    s_kwargs = (
        ylabel=raw"$\frac{S}{\log_2 N}$",
        markershape=:utriangle,
        label="entropy",
        color=:green,
    )
    pF = plot(maxdims, fidelities; f_kwargs..., kwargs...)
    # pS = plot(maxdims, entropies; s_kwargs..., kwargs...)
    return pF
end

