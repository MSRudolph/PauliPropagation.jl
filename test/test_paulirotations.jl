using Random
using Test

@testset "Test Pauli gates apply function" begin
  """Test apply function for Pauli gates."""
  nq = 1
  th = randn()
  # single qubit X gate
  gate = PauliRotation([:X], [1])
  # apply to Z
  pstr = PauliString(nq, :Z, 1)
  @test all(apply(gate, pstr.term, th) .≈ (0x03, cos(th), 0x02, sin(th)))
  # apply to Y
  pstr = PauliString(nq, :Y, 1)
  # YT: the ordering of Paulis after applying a gate is not unique.
  @test all(apply(gate, pstr.term, th) .≈ (0x02, cos(th), 0x03, -sin(th)))
  # apply to X
  pstr = PauliString(nq, :X, 1)
  @test all(apply(gate, pstr.term, th) .≈ (0x01, 1.0))

  # two-qubit Y gate
  gate = PauliRotation([:Y, :Y], [1, 2])
  # apply to Z
  pstr = PauliString(nq, :Z, 2)
  @test all(apply(gate, pstr.term, th) .≈ (0x0c, cos(th), 0x06, -sin(th)))
  # apply to Y
  pstr = PauliString(nq, :Y, 2)
  @test all(apply(gate, pstr.term, th) .≈ (0x08, 1.0))
  # apply to X
  pstr = PauliString(nq, :X, 2)
  @test all(apply(gate, pstr.term, th) .≈ (0x04, cos(th), 0x0e, sin(th)))

  # single qubit Z gate on two qubits
  nq = 2
  gate = PauliRotation([:Z], [1])
  # apply to Z
  pstr = PauliString(nq, :Z, 1)
  @test all(apply(gate, pstr.term, th) .≈ (0x03, 1.0))
  # apply to Y
  pstr = PauliString(nq, :Y, 1)
  @test all(apply(gate, pstr.term, th) .≈ (0x02, cos(th), 0x01, sin(th)))
  # apply to X
  pstr = PauliString(nq, :X, 1)
  @test all(apply(gate, pstr.term, th) .≈ (0x01, cos(th), 0x02, -sin(th)))
end

@testset "Test fastgate" begin
  nq = 2
  th = randn()
  gate = PauliRotation([:X, :Z], [1, 2])
  fastgate = tofastgates(gate, nq)
  pstr = PauliString(nq, :Z, 2)
  @test apply(gate, pstr.term, th) == apply(fastgate, pstr.term, th)
end
