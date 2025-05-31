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
                   pairs

function PauliPropagation.applytoall!(gate::PauliRotation, theta, psum::Dictionary{K,V}, aux_psum::Dictionary{K,V}) where {K,V}
    gate = _tomaskedpaulirotation(gate, K)
    cos_val = cos(theta)
    sin_val = sin(theta)
    to_update = Vector{K}()
    new_terms = Vector{Tuple{K,V}}()
    for (pstr, coeff) in pairs(psum)
        if commutes(gate, pstr)
            continue
        end
        push!(to_update, pstr)
        new_pstr, sign = getnewpaulistring(gate, pstr)
        push!(new_terms, (new_pstr, coeff * sin_val * sign))
    end
    for pstr in to_update
        psum[pstr] = psum[pstr] * cos_val
    end
    for (new_pstr, new_coeff) in new_terms
        if haskey(aux_psum, new_pstr)
            aux_psum[new_pstr] += new_coeff
        else
            insert!(aux_psum, new_pstr, new_coeff)
        end
    end
end
end
