module PauliPropagation

using Base.Threads

include("datatypes.jl")
export PathProperties, NumericPathProperties


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
    getelement,
    setelement!,
    show,
    containsXorY,
    containsYorZ

include("truncations.jl")

include("Propagation/Propagation.jl")
export mergingbfs, applygatetoall!, applygatetoone!

include("stateoverlap.jl")
export
    overlapbyorthogonality,
    overlapwithzero,
    overlapwithplus,
    orthogonaltozero,
    orthogonaltoplus,
    filterdict,
    zerofilter,
    evaluateagainstdict,
    getnumcoeff

include("surrogate.jl")
export operatortopathdict, PauliGateNode, gettraceevalorder, expectation, resetnodes

end
