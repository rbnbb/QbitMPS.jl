int2binary(x, n_bits) = Tuple(map(bit -> 1 & (x >> bit), n_bits-1:-1:0))

"""
    list_fock_states(n_qubits::Int)

Enumerate all Fock states given n qbits. Order them as eg. for 2 qbits: |00>, |01>, |10>, |11>.
"""
function list_fock_states(n_qubits::Int)
    # we can obtain the desired ordering using binary representation
    states = [int2binary(x, n_qubits) for x in 0:(2^n_qubits-1)]
    ##now put them in good order given by is
    #ordering =
    #    map(x -> parse(Int, string(tags(x))[15:end-1]), is)::Tuple{Vararg{Int,n_qubits}}
    #function swapped!(line::AbstractVector, i::Int, j::Int)
    #    line[i], line[j] = line[j], line[i]
    #end
    #for (j, ord) in enumerate(ordering)
    #    if j != ord  # swap their order
    #        println("Rearrangining...")
    #        map(line -> swapped!(line, j, ord), states)
    #    end
    #end
    #states = Tuple(states)  # because Tuples are immutable
    return states
end

"""
Return the amplitude statevector of MPS.

Useful for checking result against exact diagonalisation.
"""
function mps2statevector(ψ::Union{Vector{ITensor},MPS})::Vector{ComplexF64}
    if typeof(ψ)== MPS
        ψ = ψ.data
    end
    ITensors.set_warn_order(25)
    T = contract(ψ)
    # T = T * tags2delta(inds(T), "l=0", "l=$(N+1)")
    sv = Vector{ComplexF64}()  # statevector

    for inds in list_fock_states(length(ψ))
        inds_from_1 = CartesianIndex(Tuple(map(x -> x + 1, inds)))
        push!(sv, T[inds_from_1])
    end
    return sv
end
