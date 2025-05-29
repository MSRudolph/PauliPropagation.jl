| **Documentation**|
|:----------------:|
|[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://msrudolph.github.io/PauliPropagation.jl/stable/)[![](https://img.shields.io/badge/docs-dev-green.svg)](https://msrudolph.github.io/PauliPropagation.jl/dev/)|

# PauliPropagation.jl
`PauliPropagation.jl` is a Julia package for simulating Pauli propagation in quantum circuits and systems. It focuses on simulating the evolution of observables expressed in the Pauli basis under the action of unitary gates and non-unitary channels in a quantum circuit.

Unlike traditional simulators which simulate a circuit $\mathcal{E}$ evolving the state $\rho$ in the Schrödinger picture, Pauli propagation often adopts the Heisenberg picture, evolving an observable $O$ under $\mathcal{E}^\dagger$. This can be particularly efficient when the observables remain sparse or structured under evolution, and is useful for estimating expectation values such as $\text{Tr}\left[\rho \mathcal{E}^{\dagger}(O)\right]$, studying operator dynamics, and computing correlation functions.

Pauli propagation is related to the so-called (extended) stabilizer simulation, but is fundamentally different from, for example, tensor networks. It offers a distinct approach that can handle different regimes of quantum dynamics.

Implemented in Julia, `PauliPropagation.jl` combines high-performance computation (using features such as multiple dispatch) with an accessible and high-level interface.  

## UnitaryHACK 2025 - May 28th to June 11th
From **May 28th to June 11th** we are participating in the 5th edition of the hackathon. Solve [these GitHub Issues](https://unitaryhack.dev/projects/pauli-propagation/) and collect real money bounties &#128176;&#128176;&#128176;!

The relevant Issues are:
- For $200: https://github.com/MSRudolph/PauliPropagation.jl/issues/88
- For $150: https://github.com/MSRudolph/PauliPropagation.jl/issues/86
- For $100: https://github.com/MSRudolph/PauliPropagation.jl/issues/87
- For $50: https://github.com/MSRudolph/PauliPropagation.jl/issues/89


## Installation

> Note the current package requires `Julia 1.10+`.

The `PauliPropagation.jl` package is registered and can be installed into your environment in the following way:
```julia
using Pkg
Pkg.add("PauliPropagation")
```

### Install from GitHub
If you want to install the latest code, you can install the package directly from the Github link.
For example, if you are working with a Jupyter notebook, run
```julia
using Pkg
Pkg.add(url="https://github.com/MSRudolph/PauliPropagation.jl.git", rev="branchname")
```
where you can use the keyword `rev="branchname"` to install development versions of the package.
We don't recommend using branches other than `main` or `dev`.

### Clone the repository and install locally 
Navigate to a local directory where you want to clone this repository into and run the following in a terminal
```bash
git clone git@github.com:MSRudolph/PauliPropagation.jl.git
```
Inside this cloned repository, you can now freely import `PauliPropagation` or install it into your environment.\
Alternatively, you can push the relative path to the cloned repository to the Julia package load path called `LOAD_PATH` via
```julia
rel_path = "your/relative/path/PauliPropagation"
push!(LOAD_PATH,rel_path);
```
This may require that you have no global installation of `PauliPropagation` in your enviroment.

### A note on installing Julia 
It is recommended to install julia using `juliaup` with instructions from [here](https://github.com/JuliaLang/juliaup). Then, Julia's _long-term support_ version (currently a `1.10` version) can be installed via

```juliaup add lts```

To get started running Jupyter notebooks, start a Julia session and install the `IJulia` package.

If you are working on several projects with potentially conflicting packages, it is recommended to work with within local environments or projects.

For more details, we refer to this useful [guide](https://modernjuliaworkflows.org/writing/).

## Quick Start

You can find detailed example notebooks in the `examples` folder. We provide a brief example of how to use `PauliPropagation.jl`.

Consider simulating the dynamics of an operator $O=Z_{16}$ under the evolution of a unitary  channel $\mathcal{E}(\cdot) = U^\dagger \cdot U$ in a $n=32$ qubits system. 

```julia
using PauliPropagation

nqubits = 32

observable = PauliString(nqubits, :Z, 16) # I...IZI...I
```

Our goal is to compute

```math
\text{Tr}[U^\dagger O U \rho].
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
nlayers = 32 # l as above

circuit = tfitrottercircuit(nqubits, nlayers; topology=topology)
```

In our simulations, we can choose the circuit parameter $dt$

```julia
dt = 0.1 # time step

parameters = ones(countparameters(circuit)) * dt # all parameters
```
**Important:** The circuit and parameters are defined in the order that they would act in the Schrödinger picture. Within our `propagate()` function, the order will be _reversed_ to act on the observable. 

During the propagation via `propagate()`, we employ truncation strategies such as coefficient or weight truncations, these options can be specified as keywords. 

```julia
## the truncations
max_weight = 6 # maximum Pauli weight

min_abs_coeff = 1e-4 # minimal coefficient magnitude

## propagate through the circuit
pauli_sum = propagate(circuit, observable, parameters; max_weight, min_abs_coeff)
```
The output `pauli_sum` gives us an approximation of propagated Pauli strings

```math
U^\dagger O U \approx \sum_{\alpha} c_{\alpha} P_{\alpha}
```

Finally we can compute expectation values with an initial state such as $\rho = (|0 \rangle  \langle 0 |)^{\otimes n}$
```julia
## overlap with the initial state
overlapwithzero(pauli_sum)
# yields 0.154596728241...
```

This computation is efficient because the initial state can be written in terms of only $\mathbb{I}$ and $Z$ strings

```math
\rho = \left(\frac{\mathbb{I} + Z}{2}\right)^{\otimes n}
```

Therefore, the trace is equivalent to the sum over the coefficients of Pauli strings containing only `I` and `Z` Paulis, 

```math
\mathrm{Tr}[U^\dagger O U \rho] \approx \sum_{\alpha \in \{\mathbb{I}, Z\}\, \text{strings}} c_{\alpha}.
```

## Important Notes and Caveats
- Circuits are specified in the _Schrödinger_ picture, as if operated upon states. Behind the scenes, `propagate()` will (by default) apply the _adjoint_ circuit upon the passed `PauliSum` which is treated as the observable operator.
- Schrödinger propagation is planned but not yet supported _except_ through manually passing the _adjoint_ of the intended circuit to `propagate()`. This is often easy. For instance, with the circuit order reversed, angles in `PauliRotation` gates are negated, and `CliffordGate` are passed to `transposecliffordmap()`.
- While Pauli propagation can, in principle, be used for _extended_ stabilizer simulation, we do not currently support sub-exponential strong simulation of stabilizer states.
- Sampling quantum states is currently not supported.
- Many underlying data structures and functions can be used for other purposes involving Pauli operators.

All of the above can be addressed by writing the additional missing code due to the nice extensibility of Julia.

## Upcoming Features
This package is still work-in-progress. You will probably find certain features that you would like to have and that are currently missing.\
Here are some features that we want to implement in the future. Feel free to contribute!
- **Multi-threading and improved scalability**. Currently, PauliPropagation.jl works uses a single CPU thread and may run your hardware out of memory. Future versions should be even faster and with options to trade-off computational runtime and memory requirements. 
- **Easier Schrödinger picture propagation**. Currently, the default is Heisenberg and there is no easy way to transpose the gates.
- **A fast and flexible Surrogate version**. Currently, we provide a version of the Pauli propagation Surrogate that is _good_ and _works_, at least for Pauli gates and Clifford gates. Stay tuned for a whole lot more.

## How to contribute
We have a Slack channel `#pauli-propagation` in the [Julia Slack](https://join.slack.com/t/julialang/shared_invite/zt-2zljxdwnl-kSXbwuwFHeERyxSD3iFJdQ).

If something bothers you or you want to propose an enhancement, please open an [Issue](https://github.com/MSRudolph/PauliPropagation.jl/issues) describing everything in detail.

For a concrete change of code, please fork this GitHub repository and submit a [Pull Request](https://github.com/MSRudolph/PauliPropagation.jl/pulls).

Otherwise, feel free to reach out to the developers!

## Authors

The main developer of this package is [Manuel S. Rudolph](https://github.com/MSRudolph) in the Quantum Information and Computation Laboratory of Prof. Zoë Holmes at EPFL, Switzerland.
Contact Manuel via manuel.rudolph@epfl.ch.

Further contributors to this package include [Yanting Teng](https://github.com/teng10), [Tyson Jones](https://github.com/TysonRayJones), and [Su Yeon Chang](https://github.com/sychang42).
This package is the derivative of ongoing work at the Quantum Information and Computation lab at EPFL, supervised by Prof. Zoë Holmes.

For more specific code issues, bug fixes, etc. please open a [GitHub issue](https://github.com/MSRudolph/PauliPropagation.jl/issues).

If you are publishing research using `PauliPropagation.jl`, please cite this library and our upcoming paper presenting it (coming soon(ish)).

## Related publications
Some of the developers of this package are co-authors in the following papers using Pauli propagation and (at least parts of) this code:
- [Classical simulations of noisy variational quantum circuits](https://arxiv.org/abs/2306.05400)
- [Classical surrogate simulation of quantum systems with LOWESA](https://arxiv.org/abs/2308.09109)
- [Quantum Convolutional Neural Networks are (Effectively) Classically Simulable](https://arxiv.org/abs/2408.12739)
- [Classically estimating observables of noiseless quantum circuits](https://arxiv.org/abs/2409.01706)
- [Efficient quantum-enhanced classical simulation for patches of quantum landscapes](https://arxiv.org/abs/2411.19896)
- [Simulating quantum circuits with arbitrary local noise using Pauli Propagation](https://arxiv.org/abs/2501.13101)
  
And more are coming up.
