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
    add!,
    set!,
    empty!,
    inttosymbol,
    symboltoint,
    inttostring,
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


include("Gates/Gates.jl")
export
    Gate,
    ParametrizedGate,
    StaticGate,
    PauliRotation,
    MaskedPauliRotation,
    CliffordGate,
    clifford_map,
    reset_clifford_map!,
    createcliffordmap,
    ParametrizedNoiseChannel,
    PauliNoise,
    DepolarizingNoise,
    DephasingNoise,
    PauliXNoise,
    PauliYNoise,
    PauliZNoise,
    AmplitudeDampingNoise,
    FrozenGate,
    freeze,
    TGate

include("circuits.jl")
export
    countparameters,
    bricklayertopology,
    staircasetopology,
    get2dtopology,
    get2dstaircasetopology,
    hardwareefficientcircuit,
    efficientsu2circuit,
    tfitrottercircuit,
    tiltedtfitrottercircuit,
    heisenbergtrottercircuit,
    su4ansatz,
    qcnnansatz,
    appendSU4!,
    ibmeagletopology

include("PathProperties/PathProperties.jl")
export
    PathProperties,
    PauliFreqTracker,
    wrapcoefficients

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
