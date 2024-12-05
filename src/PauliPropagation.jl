module PauliPropagation

using Base.Threads

include("./PauliAlgebra/PauliAlgebra.jl")
export
    PauliSum,
    PauliString,
    paulis,
    coefficients,
    paulitype,
    coefftype,
    numcoefftype,
    getcoeff,
    topaulistrings,
    add,
    add!,
    subtract,
    subtract!,
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


include("pathproperties.jl")
export
    PathProperties,
    NumericPathProperties,
    wrapcoefficients

include("Gates/Gates.jl")
export
    Gate,
    ParametrizedGate,
    StaticGate,
    PauliGate,
    FastPauliGate,
    tofastgates,
    tofastgates!,
    apply,
    applynoncummuting,
    CliffordGate,
    clifford_map,
    reset_clifford_map!,
    createcliffordmap,
    applywithmap,
    ParametrizedNoiseChannel,
    PauliNoise,
    DepolarizingNoise,
    DephasingNoise,
    PauliXNoise,
    PauliYNoise,
    PauliZNoise,
    AmplitudeDampingNoise,
    FrozenGate,
    freeze


include("circuits.jl")
export
    countparameters,
    bricklayertopology,
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

include("truncations.jl")
export
    truncatedampingcoeff

include("Propagation/Propagation.jl")
export
    propagate,
    propagate!,
    mergingapply!,
    applygatetoall!,
    applygatetoone!,
    mergeandclear!,
    merge

include("stateoverlap.jl")
export
    overlapbyorthogonality,
    overlapwithzero,
    overlapwithplus,
    orthogonaltozero,
    orthogonaltoplus,
    overlapwithpaulisum,
    overlapwithmaxmixed,
    filter,
    filter!,
    zerofilter,
    zerofilter!,
    plusfilter,
    plusfilter!,
    evaluateagainstdict,
    getnumcoeff

include("numericalcertificates.jl")
export
    estimateaverageerror,
    estimateaverageerror!,
    montecarlopropagation,
    mcapply

include("Surrogate/Surrogate.jl")
export
    NodePathProperties,
    EvalEndNode,
    PauliGateNode,
    evaluate!,
    reset!

end
