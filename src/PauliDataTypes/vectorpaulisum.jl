###
##
# A file to define a Pauli sum consisting of a vector of terms and a vector of coefficients.
# Can be used for multithreaded CPU and GPU propagation.
##
###

# This type multi-threads where possible. 
using AcceleratedKernels
const AK = AcceleratedKernels


struct VectorPauliSum{TV,CV} <: AbstractPauliSum
    nqubits::Int
    terms::TV
    coeffs::CV

    function VectorPauliSum(nqubits::Int, terms::TV, coeffs::CV) where {TV,CV}
        @assert length(terms) == length(coeffs) "Length of terms and coeffs must be the same. Got $(length(terms)) and $(length(coeffs))."
        return new{TV,CV}(nqubits, terms, coeffs)
    end
end

VectorPauliSum(nqubits::Int) = VectorPauliSum(Float64, nqubits)
VectorPauliSum(::Type{CT}, nqubits::Int) where {CT} = VectorPauliSum(nqubits, getinttype(nqubits)[], CT[])

PropagationBase.storage(vpsum::VectorPauliSum) = (paulis(vpsum), coefficients(vpsum))


"""
    nqubits(vpsum::VectorPauliSum)

Get the number of qubits that the `VectorPauliSum` is defined on.
"""
nqubits(vpsum::VectorPauliSum) = vpsum.nqubits

"""
    paulis(vpsum::VectorPauliSum)

Get the vector of Pauli strings from a `VectorPauliSum`.
"""
paulis(vpsum::VectorPauliSum) = vpsum.terms

"""
    coefficients(vpsum::VectorPauliSum)
    
Get the vector of coefficients from a `VectorPauliSum`.
"""
PropagationBase.coefficients(vpsum::VectorPauliSum) = vpsum.coeffs


"""
    topaulistrings(vpsum::VectorPauliSum)

Returns the Pauli strings in a `PauliSum` and their coefficients as a list of `PauliString`.
"""
topaulistrings(vpsum::VectorPauliSum) = [PauliString(vpsum.nqubits, pauli, coeff) for (pauli, coeff) in zip(vpsum.terms, vpsum.coeffs)]


Base.similar(vpsum::VectorPauliSum) = VectorPauliSum(vpsum.nqubits, similar(vpsum.terms), similar(vpsum.coeffs))

function Base.resize!(vpsum::VectorPauliSum, n_new::Int)
    resize!(vpsum.terms, n_new)
    resize!(vpsum.coeffs, n_new)
    return vpsum
end


function Base.show(io::IO, vecpsum::VectorPauliSum)
    n_paulis = length(vecpsum)
    if n_paulis == 0
        println(io, "Empty VectorPauliSum.")
        return
    elseif n_paulis == 1
        println(io, "VectorPauliSum with 1 term:")
    else
        println(io, "VectorPauliSum with $(n_paulis) terms:")
    end

    for i in 1:length(vecpsum)
        if i > 20
            println(io, "  ...")
            break
        end
        pauli_string = inttostring(vecpsum.terms[i], vecpsum.nqubits)
        if length(pauli_string) > 20
            pauli_string = pauli_string[1:20] * "..."
        end
        println(io, vecpsum.coeffs[i], " * $(pauli_string)")
    end
end



function Base.sort!(vpsum::VectorPauliSum; by=nothing, kwargs...)
    # instead of using sortperm, we use sort!() on an index array 
    # this is to be able to sort on any properties of the terms of coeffs 

    indices = collect(1:length(vpsum))

    # default for if "by" is not provided
    byfunc = isnothing(by) ? i -> vpsum.terms[i] : by

    AK.sort!(indices; by=byfunc, kwargs...)
    vpsum.terms .= view(vpsum.terms, indices)
    vpsum.coeffs .= view(vpsum.coeffs, indices)
    return vpsum
end
