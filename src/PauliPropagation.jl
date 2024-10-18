module PauliPropagation

using Base.Threads

include("datatypes.jl")
export
    PauliSum,
    PauliString,
    add!,
    PathProperties,
    NumericPathProperties,
    wrapcoefficients


include("Gates/Gates.jl")
export
    Gate,
    PauliGate,
    FastPauliGate,
    tofastgates,
    apply,
    applynoncummuting,
    CliffordGate,
    default_clifford_map,
    reset_clifford_map!,
    applywithmap

include("circuits.jl")
export
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
