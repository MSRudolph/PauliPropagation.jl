module DictPauliPropagationExt

using Dictionaries
using Dictionaries: Dictionary
import PauliPropagation
import PauliPropagation: applytoall!, PauliRotation, _tomaskedpaulirotation, getnewpaulistring, commutes,sizehint!

using Dictionaries: Dictionary, 
                   istokenizable,
                   tokens,
                   gettokenvalue,
                   settokenvalue!,
                   haskey,
                   insert!,
                   pairs,
                   isinsertable

function PauliPropagation.applytoall!(gate::PauliRotation, theta,
        psum::Dictionary{K,V}, aux_psum::Dictionary{K,V}) where {K,V}
    gate = _tomaskedpaulirotation(gate, K)
    cos_val = cos(theta)
    sin_val = sin(theta)
    updates = Vector{Tuple{K,V}}()
    new_terms = Vector{Tuple{K,V}}()
    sizehint!(updates, length(psum))
    sizehint!(new_terms, length(psum))
    for (pstr, coeff) in pairs(psum)
        if commutes(gate, pstr)
            continue
        end
        push!(updates, (pstr, coeff * cos_val))
        new_pstr, sign = getnewpaulistring(gate, pstr)
        push!(new_terms, (new_pstr, coeff * sin_val * sign))
    end
    @inbounds for (pstr, new_coeff) in updates
        psum[pstr] = new_coeff
    end
    @inbounds for (new_pstr, new_coeff) in new_terms
        hadtoken, token = gettoken!(aux_psum, new_pstr)
        if hadtoken
            current_coeff = gettokenvalue(aux_psum, token)
            settokenvalue!(aux_psum, token, current_coeff + new_coeff)
        else
            settokenvalue!(aux_psum, token, new_coeff)
        end
    end
    return nothing
end
end

