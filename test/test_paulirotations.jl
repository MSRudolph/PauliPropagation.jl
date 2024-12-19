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
  @test propagate(gate, pstr, th) == PauliSum([PauliString(nq, :Z, 1, cos(th)), PauliString(nq, :Y, 1, sin(th))])
  # apply to Y
  pstr = PauliString(nq, :Y, 1)
  @test propagate(gate, pstr, th) == PauliSum([PauliString(nq, :Y, 1, cos(th)), PauliString(nq, :Z, 1, -sin(th))])
  # apply to X
  pstr = PauliString(nq, :X, 1)
  @test propagate(gate, pstr, th) == PauliSum(pstr)

  nq = 2
  # two-qubit Y gate
  gate = PauliRotation([:Y, :Y], [1, 2])
  # apply to Z
  pstr = PauliString(nq, :Z, 2)
  @test propagate(gate, pstr, th) == PauliSum([PauliString(nq, [:I, :Z], [1, 2], cos(th)), PauliString(nq, [:Y, :X], [1, 2], -sin(th))])
  # apply to Y
  pstr = PauliString(nq, :Y, 2)
  @test propagate(gate, pstr, th) == PauliSum(pstr)
  # apply to X
  pstr = PauliString(nq, :X, 2)
  @test propagate(gate, pstr, th) == PauliSum([PauliString(nq, [:I, :X], [1, 2], cos(th)), PauliString(nq, [:Y, :Z], [1, 2], sin(th))])

  # single qubit Z gate on three-qubit Pauli sums
  nq = 3
  gate = PauliRotation([:Z], [2])
  # apply to Z
  psum = PauliSum(PauliString(nq, :Z, 2))
  @test propagate(gate, psum, th) == psum
  # apply to Y
  psum = PauliSum(PauliString(nq, :Y, 2))
  @test propagate(gate, psum, th) == PauliSum([PauliString(nq, :Y, 2, cos(th)), PauliString(nq, :X, 2, sin(th))])
  # apply to X
  psum = PauliSum(PauliString(nq, :X, 2))
  @test propagate(gate, psum, th) == PauliSum([PauliString(nq, :X, 2, cos(th)), PauliString(nq, :Y, 2, -sin(th))])
end