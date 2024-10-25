function getpaulielement(oper::AbstractArray{T}, index::Int) where {T}
    return oper[index]
end

function getpaulielement(oper::Integer, index::Integer)
    return getbitpaulielement(oper, index)
end

function setpaulielement!(oper::AbstractArray{T}, index, element::T) where {T}
    oper[index] = element
    return oper
end

function setpaulielement!(oper::AbstractArray{T}, index, element::SinglePauliType) where {T}
    return setpaulielement!(oper, index, inttosymbol(element))
end

function setpaulielement!(oper::PauliStringType, index, element::Symbol)
    return setpaulielement!(oper, index, symboltoint(element))
end

function setpaulielement!(oper::PauliStringType, index, element::SinglePauliType)
    return setbitpaulielement!(oper, index, element)
end
