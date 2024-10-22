
function countweight(pauli_op::Vector{Symbol}, args...; kwargs...)
    return sum(op in (:X, :Y, :Z) for op in pauli_op)
end

function countweight(oper::Integer; kwargs...)
    return countbitweight(oper; kwargs...)
end


function countxy(oper::Integer; kwargs...)
    return countbitxy(oper; kwargs...)
end

function countyz(oper::Integer; kwargs...)
    return countbityz(oper; kwargs...)
end

containsXorY(symbs::AbstractArray{Symbol}, args...) = :X in symbs || :Y in symbs
containsXorY(int::Integer, args...) = countxy(int) > 0
containsYorZ(int::Integer, args...) = countyz(int) > 0


### All the commutation check functions

# TODO: Should we even support operations on non-integer operators?
function commutes(gate_generator::AbstractArray{T}, pauli_op::AbstractArray{T})::Bool where {T}
    return sum(!commutes(o1, o2) for (o1, o2) in zip(gate_generator, pauli_op)) % 2 == 0
end

function commutes(sym1::Symbol, sym2::Symbol)::Bool
    if sym1 == :I || sym2 == :I
        return true
    else
        return sym1 == sym2
    end
end

function commutes(sym1::Symbol, pauli_ind::Integer)::Bool
    return commutes(sym1, inttosymbol(pauli_ind))
end

function commutes(oper1::Integer, oper2::Integer)
    return bitcommutes(oper1, oper2)
end

function commutes(oper1::Dict{T, Float64}, oper2::Dict{T, Float64}, nq::Integer) where {T <: Integer} 
    comm = commutator(oper1, oper2, nq) 
    return isempty(comm)
end 

function commutator(oper1::T, oper2::T) where {T <: Integer}
    new_oper = zero(T)

    if commutes(oper1, oper2) 
        total_sign = 0. 
    else 
        total_sign = -1

        pw4 = 1 

        while oper1 > 0 || oper2 > 0
            sign, new_pauli_op = pauliprod(oper1&3, oper2&3, 1)

            new_oper += new_pauli_op*pw4
            total_sign *= sign
            
            pw4 = pw4 << 2
            oper1 = oper1 >> 2
            oper2 = oper2 >>2 
        end 
    end 
    return 2*total_sign, new_oper
end

function commutator(oper1::T, oper2::T, nq::Integer) where {T <: Integer}
    new_oper = zero(T)

    if commutes(oper1, oper2) 
        total_sign = 0. 
    else 
        total_sign, new_oper = pauliprod(oper1, oper2, nq)
    end 
    return 2*total_sign, new_oper
end

function commutator(oper1::Union{T, Dict{T, Float64}}, oper2::Union{T, Dict{T, Float64}}, nq::Integer) where {T <: Integer}
    oper1 = isa(oper1, Dict) ? oper1 : Dict(oper1 => 1.0)
    oper2 = isa(oper1, Dict) ? oper2 : Dict(oper1 => 1.0)
    
    new_oper = Dict{typeof(first(keys(oper1))), typeof(first(values(oper1)))}()

    for (pw1, coeff1) in oper1, (pw2, coeff2) in oper2 
        if !commutes(pw1, pw2) 
            sign, pauli = commutator(pw1, pw2, nq)

            if pauli ∉ keys(new_oper) 
                new_oper[pauli] = sign*coeff1*coeff2
            else 
                new_oper[pauli] = new_oper[pauli] + sign*coeff1*coeff2
            end 
        end
    end

    # Get rid of the pauli strings with zero coeffs
    new_oper = Dict(k=>v for (k,v) in new_oper if abs(v) != 0.) 

    return new_oper
end 

function commutator(oper1::PauliSum, oper2::PauliSum) 
    new_oper = commutator(oper1.op_dict, oper2.op_dict, oper1.nqubits)
    return PauliSum(oper1.nq, new_oper)
end 

function pauliprod(op1::Integer, op2::Integer, nq::Int)
    return pauliprod(op1, op2, 1:nq)
end

function pauliprod(op1::Symbol, op2::Integer) # assume that just one qubit is involved, TODO: This must change to one bit operation on all qubits
    return pauliprod(symboltoint(op1), op2, 1:1)
end

function singlepauliprod(op1::Integer, op2::Integer)
    return pauliprod(getelement(op1, 1), getelement(op2, 1), 1:1)
end

function pauliprod(op1::Integer, op2::Integer, changed_indices) # TODO: find a way to circumvent changed_indices being passed
    op3 = bitpauliprod(op1, op2)
    sign = calculatesign(op1, op2, op3, changed_indices)
    return sign, op3
end

function calculatesign(op1::Integer, op2::Integer, op3::Integer, changed_indices)
    sign = 1
    for qind in changed_indices
        sign *= generalizedlevicivita(
            getelement(op1, qind),
            getelement(op2, qind),
            getelement(op3, qind)
        )
    end
    return sign
end

const generalized_levicivita_matrix = cat(
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], # first arg is I
    [0 1im 0 0; 1 0 0 0; 0 0 0 1im; 0 0 -1im 0], # first arg is X
    [0 0 1im 0; 0 0 0 -1im; 1 0 0 0; 0 1im 0 0], # first arg is Y
    [0 0 0 1im; 0 0 1im 0; 0 -1im 0 0; 1 0 0 0]; # first arg is Z
    dims=3)

function generalizedlevicivita(n1::Integer, n2::Integer, n3::Integer)
    # acts like levicivita but yields the correct sign for products with I or P^2
    return generalized_levicivita_matrix[n1+1, n2+1, n3+1]
end