# Introduction

`PauliPropagation.jl` is a Julia package for Pauli propagation simulation of quantum circuits and quantum systems. It focuses on simulating the evolution of observables expressed in the Pauli basis, under both noiseless and noisy quantum dynamics.

Unlike traditional simulation approaches that evolve quantum states (Schrödinger picture), Pauli propagation often works in the Heisenberg picture, evolving observables like $\mathcal{E}^\dagger(\hat{O})$ rather than states $\mathcal{E}(\rho)$. This can be particularly efficient when the observables remain sparse or structured under evolution, and is useful for studying operator dynamics, estimating expectation values, and computing correlation functions.

Pauli propagation sits alongside other classical simulation techniques, such as stabilizer simulation and tensor networks, but offers a distinct approach that can handle different regimes of quantum complexity.

Implemented in Julia, `PauliPropagation.jl` combines high-performance computation (using features such as multiple dispatch) with an accessible design similar to Python, making it suitable both for large-scale quantum simulations.