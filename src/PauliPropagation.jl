module PauliPropagation

using Base.Threads

include("datatypes.jl")
export
    PauliSum,
    PauliString,
    getcoeff,
    getpaulistrings,
    add,
    add!,
    subtract,
    subtract!,
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
    DepolarizingNoise,
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

include("./PauliAlgebra/PauliAlgebra.jl")
export
    inttosymbol,
    symboltoint,
    inttostring,
    getpaulielement,
    setpaulielement!,
    show,
    countweight,
    countxy,
    countyz,
    containsXorY,
    containsYorZ,
    pauliprod,
    commutes,
    commutator

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

include("surrogate.jl")
export
    NodePathProperties,
    EvalEndNode,
    PauliGateNode,
    gettraceevalorder,
    expectation,
    resetnodes

end
