# Introduction

`PauliPropagation.jl` is a Julia package for Pauli propagation simulation of quantum circuits and quantum systems. It focuses on simulating the evolution of observables expressed in the Pauli basis, under e.g. the dynamics of quantum circuits.

Unlike traditional simulation approaches that evolve quantum states (Schrödinger picture), Pauli propagation often works in the Heisenberg picture, evolving observables like $\mathcal{E}^\dagger(\hat{O})$ rather than states $\mathcal{E}(\rho)$. This can be particularly efficient when the observables remain sparse or structured under evolution, and is useful for studying operator dynamics, estimating expectation values, and computing correlation functions.

Pauli propagation sits alongside other classical simulation techniques, such as stabilizer simulation and tensor networks, but offers a distinct approach that can handle different regimes of quantum complexity.

Implemented in Julia, `PauliPropagation.jl` combines high-performance computation (using features such as multiple dispatch) with an accessible design similar to Python.


## Quick Start

You can find detailed example notebooks in the `examples` folder. We provide a brief example of how to use `PauliPropagation.jl`.

Consider simulating the dynamics of an operator $O=Z_{16}$ under the evolution of a unitary  channel $\mathcal{E}(\cdot) = U^\dagger \cdot U$ in a $n=32$ qubits system. 

```julia
using PauliPropagation

nqubits = 32 # number of qubits

observable = PauliString(nqubits, :Z, 16) # observable I...IZI...I
```

Our goal is to compute

```math
\mathrm{Tr}[U^\dagger O U \rho].
```

A simple unitary $U$ is the brickwork circuit, composed of two qubit gates alternating neighbouring sites. We define the circuit connectivity by 

```julia
topology = bricklayertopology(nqubits; periodic=true)
```

where `periodic` specifies the boundary condition of the gates. The library has built-in circuits with e.g. a circuit containing alternating RX and RZZ Pauli gates on the topology. This can be defined by Trotterization of a transverse field Ising Hamiltonian with $l$ steps

```math
U = \prod_{a=1}^{l} \prod_{j=1}^n e^{-i dt   X_j} e^{-i dt Z_j Z_{j+1}}.
```

```julia
nlayers = 32 ## the number of layers $l$

circuit = tfitrottercircuit(nqubits, nlayers; topology=topology)
```

In our simulations, we can choose the circuit parameter $dt$

```julia
dt = 0.1 # time step

nparams = countparameters(circuit) # number of parameters

parameters = ones(nparams) * dt # total parameters
```

We can now compute the evolution using the `propagate` function. During the propagation, we employ truncation strategies such as coefficient or weight truncations, these options can be specified as keywords. 

```julia
## the truncations
max_weight = 6 # maximum Pauli weight

min_abs_coeff = 1e-4 # minimal coefficient magnitude

## propagate through the circuit
pauli_sum = propagate(circuit, observable, parameters; max_weight=max_weight, min_abs_coeff=min_abs_coeff)
```
The output `pauli_sum` gives us an approximation of propagated Paulis 

```math
U^\dagger O U \approx \sum_{\alpha} c_{\alpha} P_{\alpha}
```

Finally we can compute expectation values with an initial state such as $\rho=(|0 \rangle  \langle 0 |)^{\otimes n}$

```julia
## overlap with the initial state
overlapwithzero(pauli_sum)
# yields 0.154596728241...
```