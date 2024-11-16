# PauliPropagation.jl
`PauliPropagation.jl` is a Julia package for Pauli propagation simulation of quantum circuits and quantum systems.

This package is currently targeted at estimating expectation values, but many underlying data structures and functions may be used for other purposes involving Pauli operators. 

## Installation

The `PauliPropagation.jl` package is not yet registered. But you can install it in at least the following two ways:

### Install from GitHub
You can install the package using `Pkg` from the following link: `https://github.com/MSRudolph/PauliPropagation.jl.git`.\
For example, if you are working with a Jupyter notebook, run
```julia
using Pkg
Pkg.add("https://github.com/MSRudolph/PauliPropagation.jl.git")
```
where you can use the keyword `rev="branchname"` to install the version on a particular branch (although we don't recommend doing this).

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

## Examples

You can find example notebooks in the `examples` folder.

Here is a tiny working example where we approximately simulate the expectation value of a quantum circuit.
```julia
using PauliPropagation

## number of qubits
nqubits = 32

## define the observable
# here I...IZI...I
observable = PauliString(nqubits, :Z, 16)

## define the circuit
# the number of layers
nlayers = 32

# bricklayertopology is also the default if you don't provide any
topology = bricklayertopology(nqubits; periodic=true)

# a circuit containing RX and RZZ Pauli gates on the topology
# derived from the Trotterization of a transverse field Ising Hamiltonian
circuit = tfitrottercircuit(nqubits, nlayers; topology=topology)

# time step
dt = 0.1
# count the number of parameters
nparams = countparameters(circuit)
# define the parameter vector
parameters = ones(nparams) * dt

## the truncations
# maximum Pauli weight
max_weight = 6
# minimal coefficient magnitude
min_abs_coeff = 1e-4

## propagate through the circuit with our best (and currently only propagation method)
pauli_sum = mergingbfs(circuit, observable, parameters; max_weight=max_weight, min_abs_coeff=min_abs_coeff)

## overlap with the initial state
overlapwithzero(pauli_sum)
# yields 0.154596728241...
```

## Upcoming Features
This package is still work in progress. You will probably find certain features that you would like to have and that are currently missing.\
Here are some features that we want to implement in the future. Feel free to contribute!
- **A documentation website!**
- **Easier Schrödinger picture propagation**. Currently, the default is Heisenberg and there is no easy way to transpose the gates.
- **Gates with several parameters**. Currently, each parametrized gate takes exactly one parameter. Multi-parameter gates have to be decomposed.
- **A fast and flexible Surrogate version**. Currently, we provide a version of the Pauli propagation Surrogate that is _good_ and _works_, at least for Pauli gates and Clifford gates. Stay tuned for a whole lot more.

## How to contribute
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
And more are coming up.
