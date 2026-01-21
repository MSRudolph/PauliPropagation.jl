
# TODO: Once this is not typed, the old propagate() needs to be typed to ::PauliSum
function propagate(circuit, term_sum::AbstractTermSum, parameters=nothing; kwargs...)
    return propagate!(circuit, deepcopy(term_sum), parameters; kwargs...)
end


function propagate!(circuit, term_sum::TS, parameters=nothing; kwargs...) where TS<:AbstractTermSum
    prop_cache = propagate!(circuit, PropagationCache(term_sum), parameters; kwargs...)
    # we expect that there exists a constructor TS(prop_cache) for the back conversion
    try
        return TS(prop_cache)
    catch e
        throw(ErrorException("Could not convert the propagation cache of type $(typeof(prop_cache)) back to the original TermSum type $(TS).
         Make sure that a constructor $(TS)(::$(typeof(prop_cache))) is defined. Original error:\n $e"))
    end
end

function propagate!(circuit, prop_cache::AbstractPropagationCache, parameters=nothing; kwargs...)

    # TODO: The two methods should be one
    # if circ is actually a single gate, promote it to a list [gate]
    # similarly the theta if it is a single number
    circuit, parameters = _promotecircandthetas(circuit, parameters)

    # if thetas is nothing, the circuit must contain only StaticGates
    # also check if the length of thetas equals the number of parametrized gates
    _checkcircandthetas(circuit, parameters)

    # TODO: allow forward propagation 
    # For now we only assume backward propagation
    circuit = reverse(circuit)
    parameters = reverse(parameters)

    # A useful iteration tool
    parameter_iterator = Iterators.Stateful(parameters)

    for gate in circuit
        if isa(gate, ParametrizedGate)
            param = popfirst!(parameter_iterator)
            prop_cache = applymergetruncate!(gate, prop_cache, param; kwargs...)
        else
            prop_cache = applymergetruncate!(gate, prop_cache; kwargs...)
        end
    end

    return prop_cache
end


function applymergetruncate!(gate, prop_cache::AbstractPropagationCache, args...; kwargs...)

    # args is usually expected to be empty or contain a parameter for the gate
    prop_cache = applytoall!(gate, prop_cache, args...; kwargs...)

    # usually this merges from some auxillary term sum into the main term sum
    # for vector-based caches, it deduplicates within the main term sum
    prop_cache = merge!(prop_cache; kwargs...)

    prop_cache = truncate!(prop_cache; kwargs...)

    return prop_cache
end

"""
    applytoall!(gate, theta psum, output_psum; kwargs...)

2nd-level function below `applymergetruncate!` that applies one gate to all Pauli strings in `psum`, moving results into `output_psum` by default.
After this functions, Pauli strings in remaining in `psum` and `output_psum` are merged.
This function can be overwritten for a custom gate if the lower-level functions `applyandadd!` and `apply` are not sufficient.
In particular, this function can be used to manipulate both `psum` and `output_psum` at the same time to reduce memory movement.
Note that manipulating `psum` on anything other than the current Pauli string will likely lead to errors.
See the `4-custom-gates.ipynb` for examples of how to define custom gates.
"""
function applytoall!(gate, prop_cache::AbstractPropagationCache, args...; kwargs...)
    term_sum = mainsum(prop_cache)
    aux_term_sum = auxsum(prop_cache)
    # Loop over all Pauli strings in psum and apply the gate to them.
    for (pstr, coeff) in term_sum
        # apply gate to one Pauli string, move new Pauli strings to aux_psum
        applyandadd!(gate, aux_term_sum, pstr, coeff, args...; kwargs...)
    end

    # Empty psum because everything was moved into aux_term_sum. They will later be swapped.
    # If we want to reduce unnecessary Pauli string movement, we can overload applygatetoall!()
    empty!(term_sum)

    return prop_cache
end


"""
    applyandadd!(gate, pstr, coefficient, theta, output_psum; kwargs...)

3rd-level function below `applymergetruncate!` that applies one gate to one Pauli string in `psum`, moving results into `output_psum` by default.
This function can be overwritten for a custom gate if the lower-level function `apply` is not sufficient. 
This is likely the the case if `apply` is not type-stable because it does not return a unique number of outputs. 
"""
@inline function applyandadd!(gate, output_psum, pstr, coeff, args...; kwargs...)
    # Get the (potentially new) pauli strings and their coefficients in the form of ((pstr1, coeff1), (pstr2, coeff2), ...)
    pstrs_and_coeffs = apply(gate, pstr, coeff, args...; kwargs...)

    for (new_pstr, new_coeff) in pstrs_and_coeffs
        # Itererate over the pairs of pstr and coeff
        # Store the new_pstr and coeff in the aux_psum, add to existing coeff if new_pstr already exists there
        add!(output_psum, new_pstr, new_coeff)
    end

    return output_psum
end

"""
    apply(gate::StaticGate, term, coeff; kwargs...)
    apply(gate::ParametrizedGate, term, coeff; kwargs...)

Lowest-level function that applies one gate to one term and its coefficient.
Is expected to return a tuple of (new_term, new_coeff) pairs.
This function must be overloaded for each custom gate type.
Common mistakes are to return a single pair instead of a tuple of pairs, 
such as `(new_term, new_coeff)`, instead of `((new_term, new_coeff),)`.

Example:
```julia
function apply(gate::NewStaticGate, term::PauliString, coeff::Number)
    # implementation of how Hadamard gate acts on a Pauli string
    ...
    return ((new_term1, new_coeff1), (new_term2, new_coeff2), ...)
end
function apply(gate::NewParametrizedGate, term::PauliString, coeff::Number, theta::Number)
    # implementation of how RZ gate acts on a Pauli string
    ...
    return ((new_term1, new_coeff1), (new_term2, new_coeff2), ...)
end
```
"""
apply(gate, args...; kwargs...) = _thrownotimplemented(gate, :apply)




function _promotecircandthetas(circ, thetas)
    # if users pass a gate, we assume that thetas also requires a `[]` around it
    if circ isa Gate
        circ = [circ]

        if !isnothing(thetas)
            thetas = [thetas]
        end
    end

    if isnothing(thetas)
        thetas = []
    end

    return circ, thetas
end

function _checkcircandthetas(circ, thetas)
    nparams = countparameters(circ)

    if nparams != length(thetas)
        throw(ArgumentError(
            "The number of thetas must match the number of parametrized gates in the circuit. " *
            "countparameters(circ)=$nparams, length(thetas)=$(length(thetas)).")
        )
    end

    return
end