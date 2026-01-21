
function propagate(circuit, term_sum::AbstractTermSum, parameters=nothing; kwargs...)
    return propagate!(circuit, deepcopy(term_sum), parameters; kwargs...)
end


function propagate!(circuit, term_sum::TS, params=nothing; kwargs...) where TS<:AbstractTermSum
    prop_cache = propagate!(circuit, PropagationCache(term_sum), params; kwargs...)
    # we expect that there exists a constructor TS(prop_cache) for the back conversion
    # by default implemented for some storage types
    return TS(prop_cache)
end

"""
    propagate!(circuit, prop_cache::AbstractPropagationCache, params=nothing; kwargs...)

Propagate a term sum through under the action of gates in `circuit` using a propagation cache `prop_cache`. 
This is in-place and modifies `prop_cache` along with the term sums it carries.
Parameters for the parametrized gates in `circuit` are given by `params`, and need to be passed as if the circuit was applied as written in the Schrödinger picture.
If params are not passed, the circuit must contain only non-parametrized `StaticGates`.
`kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, `applyandadd!`, and `apply`,
as well as `merge!` and `truncate!`.
Default truncation kwargs are `min_abs_coeff` and `customtruncfunc`.
"""
function propagate!(circuit, prop_cache::AbstractPropagationCache, params=nothing; kwargs...)

    # TODO: The two methods should be one
    # if circuit is actually a single gate, promote it to a list [gate]
    # similarly the params if it is a single number
    circuit, params = _promotecircandparams(circuit, params)

    # if params is nothing, the circuit must contain only StaticGates
    # also check if the length of params equals the number of parametrized gates
    _checkcircandparams(circuit, params)

    # TODO: allow forward propagation 
    # For now we only assume backward propagation
    circuit = reverse(circuit)
    params = reverse(params)

    # A useful iteration tool
    parameter_iterator = Iterators.Stateful(params)

    for gate in circuit
        if isa(gate, ParametrizedGate)
            param = popfirst!(parameter_iterator)
            applymergetruncate!(gate, prop_cache, param; kwargs...)
        else
            applymergetruncate!(gate, prop_cache; kwargs...)
        end
    end

    return prop_cache
end

"""
    applymergetruncate!(gate, prop_cache::AbstractPropagationCache; kwargs...)
    applymergetruncate!(gate, prop_cache::AbstractPropagationCache, parameter; kwargs...)

1st-level function below `propagate!` that applies one gate to all terms in the main term sum `mainsum(prop_cache)`, 
potentially using an auxiliary term sum `auxsum(prop_cache)` in the process. 
All terms are then merged and deduplicated into the main term sum.
Truncations are performed after merging.
This function can be overwritten for a custom gate if the lower-level functions `applytoall!`, `applyandadd!`, and `apply` are not sufficient.
"""
function applymergetruncate!(gate, prop_cache::AbstractPropagationCache, args...; kwargs...)

    # args is usually expected to be empty or contain a parameter for the gate
    # prop_cache is modified in place
    applytoall!(gate, prop_cache, args...; kwargs...)

    # usually this merges from some auxillary term sum into the main term sum
    # for vector-based caches, it deduplicates within the main term sum
    merge!(prop_cache; kwargs...)

    truncate!(prop_cache; kwargs...)

    return
end

"""
    applytoall!(gate, prop_cache::AbstractPropagationCache; kwargs...)
    applytoall!(gate, prop_cache::AbstractPropagationCache, parameter; kwargs...)

1st-level function below `propagate!` that applies one gate to all terms in the main term sum `term_sum = mainsum(prop_cache)`, 
potentially using an auxiliary term sum `aux_term_sum = auxsum(prop_cache)` in the process. 
After this functions, all terms remaining in `term_sum` and `aux_term_sum` are merged.
This function can be overwritten for a custom gate if the lower-level functions `applyandadd!` and `apply` are not sufficient.
In particular, this function can be used to manipulate both `term_sum` and `aux_term_sum` at the same time to reduce memory movement.
Note that manipulating `term_sum` on anything other than the current term will likely lead to errors.
"""
function applytoall!(gate, prop_cache::AbstractPropagationCache, args...; kwargs...)
    term_sum = mainsum(prop_cache)
    aux_term_sum = auxsum(prop_cache)
    # Loop over all terms in term_sum and apply the gate to them.
    for (term, coeff) in term_sum
        # apply gate to one term, move new terms to aux_term_sum
        applyandadd!(gate, aux_term_sum, term, coeff, args...; kwargs...)
    end

    # Empty term_sum because everything was moved into aux_term_sum. They will later be swapped.
    # If we want to reduce unnecessary Pauli string movement, we can overload applygatetoall!()
    empty!(term_sum)

    # by default we can already swap the term sums here
    # merge!(prop_cache) will then likely not do anything
    swapsums!(prop_cache)

    return
end


"""
    applyandadd!(gate::StaticGate, output_term_sum, term, coeff; kwargs...)
    applyandadd!(gate::ParametrizedGate, output_term_sum, term, coeff, param; kwargs...)

3rd-level function below `propagate!` that applies a gate to one term in `term_sum`, moving results into `output_term_sum` by default,
which is modified in place.
This function can be overwritten for a custom gate if the lower-level function `apply` is not sufficient. 
This is likely the the case if `apply` is not type-stable because it does not return a unique number of outputs. 
"""
@inline function applyandadd!(gate, output_term_sum, term, coeff, args...; kwargs...)
    # Get the (potentially new) terms and their coefficients in the form of ((term1, coeff1), (term2, coeff2), ...)
    terms_and_coeffs = apply(gate, term, coeff, args...; kwargs...)

    for (new_term, new_coeff) in terms_and_coeffs
        # Itererate over the pairs of term and coeff
        # Store the new_term and coeff in the aux_term_sum, add to existing coeff if new_term already exists there
        add!(output_term_sum, new_term, new_coeff)
    end

    return
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
function apply(gate::NewStaticGate, term, coeff)
    # how some gate acts in the term and updates coeff
    ...
    return ((new_term1, new_coeff1), (new_term2, new_coeff2), ...)
end
function apply(gate::NewParametrizedGate, term, coeff, param)
    # how some gate acts in the term and updates coeff
    # a parameter `param` will be automatically passed if NewParametrizedGate <: ParametrizedGate
    ...
    return ((new_term1, new_coeff1), (new_term2, new_coeff2), ...)
end
```
"""
@inline apply(gate, args...; kwargs...) = _thrownotimplemented(gate, :apply)




function _promotecircandparams(circ, params)
    # if users pass a gate, we assume that thetas also requires a `[]` around it
    if circ isa Gate
        circ = [circ]

        if !isnothing(params)
            params = [params]
        end
    end

    if isnothing(params)
        params = []
    end

    return circ, params
end

function _checkcircandparams(circ, params)
    nparams = countparameters(circ)

    if nparams != length(params)
        throw(ArgumentError(
            "The number of parameters must match the number of parametrized gates in the circuit. " *
            "countparameters(circ)=$nparams, length(params)=$(length(params)).")
        )
    end

    return
end