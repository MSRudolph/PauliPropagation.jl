### pathweighttracker.jl
##
# This file contains the PathWeightTracker type and specializations for propagation.
# Whenever a PauliSum with PathWeightTracker coefficients is propagated through an AmplitudeDampingNoise or PauliNoise gate,
# and when it does not act on the identity, the path weight is incremented by one.
# This can be used within custom truncation functions, for example.
# Wrap coefficients into PathWeightTracker by using the `wrapcoefficients()` function.
##
###
struct PathWeightTracker <: PathProperties
    coeff::Float64
    path_weight::Int
end

PathWeightTracker(coeff) = PathWeightTracker(coeff, 0)

function PauliPropagation.applytoall!(gate::AmplitudeDampingNoise, gamma, psum::PauliSum{TT,PathWeightTracker}, aux_psum::PauliSum{TT,PathWeightTracker}, args...; kwargs...) where TT

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum

        pauli = getpauli(pstr, gate.qind)

        if pauli == 0
            # Pauli is I, so the gate does not do anything
            continue
        end

        coeff = PathWeightTracker(coeff.coeff, coeff.path_weight + 1)

        if pauli == 1 || pauli == 2
            # Pauli is X or Y, so the gate will give a sqrt(1-gamma) prefactor
            new_coeff = sqrt(1 - gamma) * coeff
            # set the coefficient of the Pauli string in the psum to the new coefficient
            set!(psum, pstr, new_coeff)
        else
            # Pauli is Z, so the gate will split the Pauli string 

            # else we know the gate will split th Pauli string into two
            new_pstr = setpauli(pstr, 0, gate.qind)
            coeff1 = (1 - gamma) * coeff
            coeff2 = gamma * coeff

            # set the coefficient of the original Pauli string
            set!(psum, pstr, coeff1)

            # add the coefficient of the new Pauli string in the aux_psum
            add!(aux_psum, new_pstr, coeff2)

        end
    end

    return
end

function PauliPropagation.applytoall!(gate::PauliNoise, p, psum::PauliSum{TT,PathWeightTracker}, aux_psum::PauliSum{TT,PathWeightTracker}; kwargs...) where TT<:PauliStringType

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        if getpauli(pstr, gate.qind) == 0
            # Pauli is I, so the gate does not do anything
            continue
        end

        # apply the Pauli noise, which will reduce the coefficient
        new_coeff = PathWeightTracker(coeff.coeff * (1 - p), coeff.path_weight + 1)

        # set the coefficient of the Pauli string in the psum to the new coefficient
        set!(psum, pstr, new_coeff)
    end

    return
end