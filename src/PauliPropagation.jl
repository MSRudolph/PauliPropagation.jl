module PauliPropagation

using Base.Threads
using LinearAlgebra

include("./PauliAlgebra/PauliAlgebra.jl")
export
    PauliStringType,
    PauliType,
    PauliSum,
    PauliString,
    paulis,
    coefficients,
    norm,
    paulitype,
    coefftype,
    numcoefftype,
    getcoeff,
    topaulistrings,
    mult!,
    add!,
    set!,
    mult!,
    empty!,
    similar,
    identitypauli,
    identitylike,
    inttosymbol,
    symboltoint,
    inttostring,
    ispauli,
    getpauli,
    setpauli,
    show,
    countweight,
    countxy,
    countyz,
    containsXorY,
    containsYorZ,
    pauliprod,
    commutes,
    commutator,
    getinttype

include("PauliTransferMatrix/PauliTransferMatrix.jl")
export
    calculateptm,
    totransfermap

include("Gates/Gates.jl")
export
    Gate,
    ParametrizedGate,
    StaticGate,
    PauliRotation,
    MaskedPauliRotation,
    CliffordGate,
    clifford_map,
    transposecliffordmap,
    reset_clifford_map!,
    createcliffordmap,
    composecliffordmaps,
    ParametrizedNoiseChannel,
    PauliNoise,
    DepolarizingNoise,
    DephasingNoise,
    AmplitudeDampingNoise,
    PauliXNoise,
    PauliYNoise,
    PauliZNoise,
    FrozenGate,
    freeze,
    TGate,
    TransferMapGate,
    tomatrix

include("Circuits/Circuits.jl")
export
    countparameters,
    getparameterindices,
    bricklayertopology,
    staircasetopology,
    rectangletopology,
    staircasetopology2d,
    ibmeagletopology,
    hardwareefficientcircuit,
    efficientsu2circuit,
    tfitrottercircuit,
    tiltedtfitrottercircuit,
    heisenbergtrottercircuit,
    su4circuit,
    qcnncircuit,
    appendSU4!,
    rxlayer!,
    rylayer!,
    rzlayer!,
    rxxlayer!,
    ryylayer!,
    rzzlayer!

include("PathProperties/PathProperties.jl")
export
    PathProperties,
    PauliFreqTracker,
    wrapcoefficients,
    unwrapcoefficients

include("truncations.jl")
export
    truncateweight,
    truncatemincoeff,
    truncatefrequency,
    truncatesins,
    truncatedampingcoeff

include("Propagation/Propagation.jl")
export
    propagate,
    propagate!,
    applymergetruncate!,
    applytoall!,
    apply,
    applyandadd!,
    mergeandempty!,
    merge

include("stateoverlap.jl")
export
    overlapbyorthogonality,
    overlapwithzero,
    overlapwithplus,
    overlapwithones,
    orthogonaltozero,
    orthogonaltoplus,
    overlapwithcomputational,
    overlapwithpaulisum,
    overlapwithmaxmixed,
    filter,
    filter!,
    zerofilter,
    zerofilter!,
    plusfilter,
    plusfilter!,
    evaluateagainstdict,
    tonumber

include("numericalcertificates.jl")
export
    estimatemse,
    estimatemse!

include("Surrogate/Surrogate.jl")
export
    NodePathProperties,
    evaluate!,
    reset!

end
