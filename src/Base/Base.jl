module PropagationBase
using LinearAlgebra
using AcceleratedKernels
const AK = AcceleratedKernels

include("./utils.jl")
export tonumber

include("./termsum.jl")
export
    AbstractTermSum,
    terms,
    coefficients,
    coeffs,
    termtype,
    coefftype,
    numcoefftype,
    getcoeff,
    nsites,
    add!,
    mult!,
    set!,
    empty!,
    similar,
    capacity

include("./propagationcache.jl")
export
    AbstractPropagationCache,
    PropagationCache,
    mainsum,
    auxsum,
    setmainsum!,
    setauxsum!,
    swapsums!,
    activesize,
    setactivesize!,
    activeterms,
    activecoeffs,
    activeauxterms,
    activeauxcoeffs,
    flags,
    indices,
    activeflags,
    activeindices,
    lastactiveindex,
    resize!

include("./gates.jl")
export
    Gate,
    StaticGate,
    ParametrizedGate,
    countparameters

include("./propagate.jl")
export propagate,
    propagate!,
    applymergetruncate!,
    applytoall!,
    applyandadd!,
    apply

include("./merge.jl")
export merge!, mergefunc

include("./truncate.jl")
export truncate!


include("./vectorbackend.jl")
export
    sortbyterm!,
    flag!,
    flagterms!,
    flagcoeffs!,
    flagstoindices!,
    permuteviaindices!,
    filterviaflags!


end