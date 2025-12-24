###
##
# A file to define a Pauli sum consisting of a vector of terms and a vector of coefficients.
# Can be used for multithreaded CPU and GPU propagation.
##
###

struct VectorPauliSum{TV,CV}
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
VectorPauliSum(pstr::PauliString) = VectorPauliSum(pstr.nqubits, [pstr.term], [pstr.coeff])
VectorPauliSum(psum::PauliSum) = VectorPauliSum(psum.nqubits, collect(paulis(psum)), collect(coefficients(psum)))

paulis(vpsum::VectorPauliSum) = vpsum.terms
coefficients(vpsum::VectorPauliSum) = vpsum.coeffs
paulitype(vpsum::VectorPauliSum) = eltype(vpsum.terms)
coefftype(vpsum::VectorPauliSum) = eltype(vpsum.coeffs)

Base.similar(vpsum::VectorPauliSum) = VectorPauliSum(vpsum.nqubits, similar(vpsum.terms), similar(vpsum.coeffs))

function Base.resize!(vpsum::VectorPauliSum, n_new::Int)
    resize!(vpsum.terms, n_new)
    resize!(vpsum.coeffs, n_new)
    return vpsum
end