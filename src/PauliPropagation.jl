module PauliPropagation

using Base.Threads
using ThreadsX

include("./PauliAlgebra/PauliAlgebra.jl")
export
    PauliSum,
    PauliString,
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
    apply,
    applynoncummuting,
    CliffordGate,
    default_clifford_map,
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
    AmplitudeDampingNoise


include("circuits.jl")
export
    countparameters,
    bricklayertopology,
    get2dtopology,
    get2dstaircasetopology,
    hardwareefficientcircuit,
    efficientsu2circuit,
    tfitrottercircuit,
    heisenbergtrottercircuit,
    su4ansatz,
    qcnnansatz,
    appendSU4!

include("truncations.jl")
export
    truncatedampingcoeff

include("Propagation/Propagation.jl")
export
    mergingbfs,
    mergingbfs!,
    applygatetoall!,
    applygatetoone!

include("stateoverlap.jl")
export
    overlapbyorthogonality,
    overlapwithzero,
    overlapwithplus,
    orthogonaltozero,
    orthogonaltoplus,
    filter,
    zerofilter,
    evaluateagainstdict,
    getnumcoeff

include("numericalcertificates.jl")
export
    estimateaverageerror,
    estimateaverageerror!,
    montecarlopropagation,
    mcapply

include("surrogate.jl")
export
    NodePathProperties,
    EvalEndNode,
    PauliGateNode,
    gettraceevalorder,
    traceevalorder,
    expectation,
    resetnodes

end
