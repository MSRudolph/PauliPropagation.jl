struct CliffordGate <: StaticGate
    symbol::Symbol
    qinds::Vector{Int}
end

# TODO: verify that these are all correct
const _default_clifford_map = Dict(
    :H => [(1, 0x00), (1, 0x03), (-1, 0x02), (1, 0x01)],
    :X => [(1, 0x00), (1, 0x01), (-1, 0x02), (-1, 0x03)],
    :Y => [(1, 0x00), (-1, 0x01), (1, 0x02), (1, 0x03)],
    :Z => [(1, 0x00), (-1, 0x01), (-1, 0x02), (1, 0x03)],
    :S => [(1, 0x00), (-1, 0x02), (1, 0x01), (1, 0x03)],
    :CNOT => [(1, 0x00), (1, 0x05), (1, 0x06), (1, 0x03), (1, 0x04), (1, 0x01), (1, 0x02), (1, 0x07), (1, 0x0b), (1, 0x0e), (-1, 0x0d), (1, 0x08), (1, 0x0f), (-1, 0x0a), (1, 0x09), (1, 0x0c)],
    :ZZpihalf => [(1, 0x00), (1, 0x0e), (-1, 0x0d), (1, 0x03), (1, 0x0b), (1, 0x05), (1, 0x06), (1, 0x08), (-1, 0x07), (1, 0x09), (1, 0x0a), (-1, 0x04), (1, 0x0c), (1, 0x02), (-1, 0x01), (1, 0x0f)],
    :SWAP => [
        (1, 0x00), (1, 0x04), (1, 0x08), (1, 0x0c), (1, 0x01), (1, 0x05), 
        (1, 0x09), (1, 0x0d), (1, 0x02), (1, 0x06), (1, 0x0a), (1, 0x0e), 
        (1, 0x03), (1, 0x07), (1, 0x0b), (1, 0x0f)
    ],
)

const default_clifford_map = deepcopy(_default_clifford_map)

function reset_clifford_map!()
    global default_clifford_map = deepcopy(_default_clifford_map)
    return
end


# Generalized function to create a gate map for any gate
function createcliffordmap(gate_relations::Dict)
  """Convert gate relations to a Clifford gate map.
  
  Args:
  gate_relations: Dict{Tuple, Tuple} - A dictionary of gate relations

  Returns:
  mapped_gate: Vector{Tuple} - A vector of tuples representing the gate map.
  """
  # check that number of qubits are less than 5 (supported by UInt8)
  if maximum([length(k) for k in keys(gate_relations)]) > 4
      throw(ArgumentError(
          "Number of qubits less than 5 is supported for UInt8 type."
      ))
  end

  gate_keys = collect(keys(gate_relations))
  order_indices = [symboltoint(collect(k)) for k in gate_keys]
  
  # Initialize arrays for reordered keys and values
  reordered_gate_vals = Vector{
      Tuple{Int, typeof(gate_keys[1]).parameters...}
  }(undef, length(gate_keys))
  for (i, idx) in enumerate(order_indices)
      reordered_gate_vals[idx + 1] = gate_relations[gate_keys[i]]
  end
  
  mapped_gate = Vector{
      Tuple{Int, typeof(symboltoint(collect(gate_keys[1])))}
  }(undef, length(gate_keys))
  for (i, v) in enumerate(reordered_gate_vals)
      mapped_gate[i] = v[1], symboltoint(collect(v[2:end]))
  end
  
  return mapped_gate
end

### Applying Clifford gates
function apply(gate::CliffordGate, operator, theta, coefficient) # TODO: write tests for this
    return apply(gate, operator, coefficient)
end

function apply(gate::CliffordGate, operator, coefficient=1.0) # TODO: write tests for this
    map_array = default_clifford_map[gate.symbol]
    return applywithmap(gate, operator, coefficient, map_array)
end

function applywithmap(gate, operator, coefficient, map_array)
    operator = copy(operator)
    qinds = gate.qinds

    lookup_op = _extractlookupop(operator, qinds)
    sign, new_op = map_array[lookup_op+1]  # +1 because Julia is 1-indexed and lookup_op is 0-indexed
    operator = _insertnewop!(operator, new_op, qinds)

    coefficient = _multiplysign!(coefficient, sign)
    return operator, coefficient
end

function _extractlookupop(operator, qinds)
    lookup_op = typeof(operator)(0)
    for ii in eachindex(qinds)
        lookup_op = setelement!(lookup_op, ii, getelement(operator, qinds[ii]))
    end
    return lookup_op
end

function _insertnewop!(operator, new_op, qinds)
    for ii in eachindex(qinds)
        operator = setelement!(operator, qinds[ii], getelement(new_op, ii))
    end
    return operator
end

function _multiplysign!(coefficient::Number, sign)
    return coefficient * sign
end

function _multiplysign!(coefficient::NumericPathProperties, sign)
    coefficient.coeff *= sign
    return coefficient
end
