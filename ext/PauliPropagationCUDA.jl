# ext/PauliPropagationCUDA.jl
module PauliPropagationCUDA

using PauliPropagation
using CUDA

# TODO: export utilities


function CUDA.cu(psum::VectorPauliSum)
    cu_paulis = cu(paulis(psum))
    cu_coeffs = cu(coefficients(psum))
    return VectorPauliSum(nqubits(psum), cu_paulis, cu_coeffs)
end

function CUDA.cu(prop_cache::VectorPauliPropagationCache)
    cu_mainsum = cu(mainsum(prop_cache))
    cu_auxsum = cu(auxsum(prop_cache))
    cu_flags = cu(prop_cache.flags)
    cu_indices = cu(prop_cache.indices)
    return VectorPauliPropagationCache(cu_mainsum, cu_auxsum, cu_flags, cu_indices, prop_cache.active_size)
end

function Base.collect(psum::VectorPauliSum)
    return VectorPauliSum(nqubits(psum), collect(paulis(psum)), collect(coefficients(psum)))
end

function Base.show(io::IO, psum::VectorPauliSum{TA,CA}) where {TA<:CUDA.CuArray,CA<:CUDA.CuArray}
    println(io, "CUDA VectorPauliSum ($(nqubits(psum)) qubits, $(length(coefficients(psum))) Paulis)")
end

function Base.show(io::IO, prop_cache::VectorPauliPropagationCache{TA,CA,FA,IA}) where {TA<:CUDA.CuArray,CA<:CUDA.CuArray,FA<:CUDA.CuArray,IA<:CUDA.CuArray}
    println(io, "CUDA VectorPauliPropagationCache with active size $(length(prop_cache)) (capacity: $(capacity(prop_cache)))")
end

function PauliPropagation.lastactiveindex(prop_cache::VectorPauliPropagationCache{TA,CA,FA,IA}) where {TA<:CUDA.CuArray,CA<:CUDA.CuArray,FA<:CUDA.CuArray,IA<:CUDA.CuArray}
    return CUDA.@allowscalar prop_cache.indices[prop_cache.active_size]
end



end