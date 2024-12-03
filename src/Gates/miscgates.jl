### miscgates.jl
##
# A file for gates that don't fit into the other categories.
##
###


## Implement a T-gate as a frozen RZ rotation gate. Easy and fast enough.
"""
    TGate(qind::Integer)

Returns a T gate acting on qubit `qind`.
This gate is implemented as a frozen `PauliGate` acting on the `qind`-th qubit with a rotation of π/4 around the Z axis.
"""
function TGate(qind)
    if !isa(qind, Integer)
        throw("TGate requires an integer site index `qind`. Got $(typeof(qind)).")
    end
    return PauliGate(:Z, qind, π / 4)
end

"""
Pretty print for `TGate` implemented as frozen Z-Pauli gate.
"""
function Base.show(io::IO, frozengate::FrozenGate{PG,Float64}) where {PG<:Union{PauliGate,FastPauliGate}}
    gate_str = frozengate.gate isa PauliGate ? "TGate" : "FastTGate"
    print(io, gate_str * "($(frozengate.gate.qinds[1]))")
end