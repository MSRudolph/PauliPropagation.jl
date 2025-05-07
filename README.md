# PauliPropagation.jl
`PauliPropagation.jl` is a Julia package for Pauli propagation simulation of quantum circuits and quantum systems. It focuses on simulating the evolution of observables expressed in the Pauli basis, under e.g. the dynamics of quantum circuits.

Unlike traditional simulation approaches that evolve quantum states (Schrödinger picture), Pauli propagation often works in the Heisenberg picture, evolving observables like $\mathcal{E}^\dagger(\hat{O})$ rather than states $\mathcal{E}(\rho)$. This can be particularly efficient when the observables remain sparse or structured under evolution, and is useful for studying operator dynamics, estimating expectation values, and computing correlation functions.

Pauli propagation sits alongside other classical simulation techniques, such as stabilizer simulation and tensor networks, but offers a distinct approach that can handle different regimes of quantum complexity.

Implemented in Julia, `PauliPropagation.jl` combines high-performance computation (using features such as multiple dispatch) with an accessible design similar to Python.
<!-- `PauliPropagation.jl` is a Julia package for Pauli propagation simulation of quantum circuits and quantum systems.

The package simulates the evolution of objects expressed in the Pauli basis under noiseless and noisy quantum circuits. Commonly, this is used for the Heisenberg picture evolution of an observable. For example, if $`\hat{O}`$ is an observable that is preferably sparse in Pauli basis and $`\mathcal{E}`$ is a quantum circuit, we simulate $`\mathcal{E}^\dagger(\hat{O})`$ instead of most quantum simulation packages simulating the Schrödinger evolution  $`\mathcal{E}(\rho)`$ of states $`\rho`$. For the case of unitary quantum circuits $`U`$, the evolved observable $`\mathcal{E}^\dagger(\hat{O})`$ is usually written like $`U^\dagger \hat{O} U`$.

Some opt-in truncations or approximations are particularly suited for estimating expectation values $`Tr[\rho \mathcal{E}^\dagger(\hat{O})]`$ of evolved observables with quantum states.  -->


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

### Clone repository and install locally 
Navigate to a local directory where you want to clone this repository into and run the following in a terminal
```bash
git clone git@github.com:MSRudolph/PauliPropagation.jl.git
```
Inside this cloned repository you can now freely import `PauliPropagation` or install it into your environment.\
Alternatively, you can push the relative path to the cloned repository to the Julia package load path called `LOAD_PATH` via
```julia
rel_path = "your/relative/path/PauliPropagation"
push!(LOAD_PATH,rel_path);
```
This may require that you have no global installation of `PauliPropagation` in your enviroment.

### A note on julia installation 
It is recommended to install julia using `juliaup`

```juliaup add 1.10```

Go to the project directory (e.g. PauliPropagation.jl). To start julia for a local environment

```julia --project=./```

More details can be found on this useful [guide](https://modernjuliaworkflows.org/writing/#installation).

## Quick Start

You can find detailed example notebooks in the `examples` folder. We provide a brief example of how to use `PauliPropagation.jl`.

Consider simulating the dynamics of an operator $O=Z_{16}$ under the evolution of a unitary  channel $\mathcal{E}(\cdot) = U^\dagger \cdot U$ in a $n=32$ qubits system. 

```julia
using PauliPropagation

nqubits = 32

observable = PauliString(nqubits, :Z, 16) # I...IZI...I
```

Our goal is to compute

$$
\text{Tr}[U^\dagger O U \rho].
$$

A simple unitary $U$ is the brickwork circuit, composed of two qubit gates alternating neighbouring sites. We define the circuit connectivity by 

```julia
topology = bricklayertopology(nqubits; periodic=true)
```

where `periodic` specifies the boundary condition of the gates. The library has built-in circuits with e.g. a circuit containing alternating RX and RZZ Pauli gates on the topology. This can be defined by Trotterization of a transverse field Ising Hamiltonian with $l$ steps

$$
U = \prod_{a=1}^{l} \prod_{j=1}^n e^{-i dt   X_j} e^{-i dt Z_j Z_{j+1}}.
$$

```julia
nlayers = 32 # l as above

circuit = tfitrottercircuit(nqubits, nlayers; topology=topology)
```

In our simulations, we can choose the circuit parameter $dt$

```julia
dt = 0.1 # time step

parameters = ones(countparameters(circuit)) * dt # all parameters
```

We can now compute the evolution using the `propagate` function. During the propagation, we employ truncation strategies such as coefficient or weight truncations, these options can be specified as keywords. 

```julia
## the truncations
max_weight = 6 # maximum Pauli weight

min_abs_coeff = 1e-4 # minimal coefficient magnitude

## propagate through the circuit
pauli_sum = propagate(circuit, observable, parameters;
max_weight=max_weight, min_abs_coeff=min_abs_coeff)
```
The output `pauli_sum` gives us an approximation of propagated Paulis 

$$
U^\dagger O U \approx \sum_{\alpha} c_{\alpha} P_{\alpha}
$$

Finally we can compute expectation values with an initial state such as $\rho = (|0 \rangle  \langle 0 |)^{\otimes n}$
```julia
## overlap with the initial state
overlapwithzero(pauli_sum)
# yields 0.154596728241...
```

<!-- This computation is efficient because the initial state can be written in terms of only $\mathbb{I}$ and $Z$ strings
$$
\rho = (\frac{\mathbb{I} + Z}{2})^{\otimes n}
$$
Therefore, the trace is equivalent to sum over coefficients of such Pauli strings 
$$
\mathrm{Tr}[U^\dagger O U \rho] \approx \sum_{\alpha \in \{\mathbb{I}, Z\} \text{strings}} c_{\alpha}.
$$ -->

## Important Notes and Caveats
All of the following points can be addressed by you writing the necessary missing code due to the nice extensibility of Julia.
- The default is the Heisenberg _backpropagation_. Schrödinger propagation may soon be natively supported. At this moment, there are options to transpose `PauliRotation` gates by multiplying their angles with `-1` and `CliffordGate`s by using `transposecliffordmap()`.
- We currently do not support the strong simulation of quantum states in non-exponential time (even for Stabilizer states). Pauli propagation could in principle be used as a backend for extended stabilizer simulation.
- Sampling quantum states is currently not supported.
- Many underlying data structures and functions can be used for other purposes involving Pauli operators.

## Upcoming Features
This package is still work-in-progress. You will probably find certain features that you would like to have and that are currently missing.\
Here are some features that we want to implement in the future. Feel free to contribute!
- **A documentation website!**
- **Easier Schrödinger picture propagation**. Currently, the default is Heisenberg and there is no easy way to transpose the gates.
- **A fast and flexible Surrogate version**. Currently, we provide a version of the Pauli propagation Surrogate that is _good_ and _works_, at least for Pauli gates and Clifford gates. Stay tuned for a whole lot more.

## How to contribute
We have a Slack channel `#pauli-propagation` in the [Julia Slack](https://join.slack.com/t/julialang/shared_invite/zt-2zljxdwnl-kSXbwuwFHeERyxSD3iFJdQ).

If something bothers you or you want to propose an enhancement, please open an [Issue](https://github.com/MSRudolph/PauliPropagation.jl/issues) describing everything in detail.

For a concrete change of code, please fork this GitHub repository and submit a [Pull Request](https://github.com/MSRudolph/PauliPropagation.jl/pulls).

Otherwise, feel free to reach out to the developers!

## Authors

The main developer of this package is Manuel S. Rudolph in the Quantum Information and Computation Laboratory of Prof. Zoë Holmes at EPFL, Switzerland.
Contact Manuel via manuel.rudolph@epfl.ch.

This package is the derivative of ongoing collaborations with Armando Angrisani and [Tyson Jones](https://github.com/TysonRayJones) at EPFL, supervised by Prof. Zoë Holmes at EPFL.

Further contributors to this package include [Yanting Teng](https://github.com/teng10) and [Su Yeon Chang](https://github.com/sychang42).

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
