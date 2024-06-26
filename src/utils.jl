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
function mps2statevector(
    ψ::Union{Vector{ITensor},MPS};
    reverse_order = false,
)::Vector{ComplexF64}
    if typeof(ψ) == MPS
        ψ = ψ.data
    end
    ITensors.set_warn_order(25)
    T = contract(ψ)
    # T = T * tags2delta(inds(T), "l=0", "l=$(N+1)")
    sv = Vector{ComplexF64}()  # statevector

    for inds in list_fock_states(length(ψ))
        if reverse_order
            product_state = Tuple(reverse(map(x -> x + 1, inds)))
        else
            product_state = Tuple(map(x -> x + 1, inds))
        end
        inds_from_1 = CartesianIndex(product_state)
        push!(sv, T[inds_from_1])
    end
    return sv
end

ppsi(psi) = [join.(list_fock_states(length(psi))) mps2statevector(psi)]

"""

SeeITensor - pretty print ITensor with the Qubit index as the third dimension,
this allows interpreting it as the notation:
[ |0> ]
"""
function seit(A::ITensor)
    sidx = inds(A; tags="Site")
    lidx = inds(A; tags="Link")
    sitenum(l) = parse(Int, match(r"[ln]=(\d)", string(tags(l))).captures[1]) # site number
    if length(lidx)>1 && sitenum(lidx[1]) > sitenum(lidx[2])  #leftmost has smaller number
        lidx = (lidx[2], lidx[1])
    end
    return permute(A, lidx..., sidx...; allow_alias=true).tensor
end

# debug helper line:
# <BS>psi = productMPS(siteinds(mpo), ["0", "1", "0", "0"]);apply!(mpo, psi); ppsi(psi)[ppsi(psi)[:,2] .!= 0, :]
