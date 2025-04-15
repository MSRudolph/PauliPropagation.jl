using Documenter, PauliPropagation

# determines doc site layout
makedocs(
    sitename="Pauli Propagation",
    pages=[
        "index.md",
        "installation.md",
        "tutorials.md",
        "API" => [
            "api/Circuits.md",
            "api/Gates.md",
            "api/PathProperties.md",
            "api/PauliAlgebra.md",
            "api/PauliTransferMatrix.md",
            "api/Propagation.md",
            "api/Surrogate.md",
            "api/NumericalCertificates.md",
            "api/StateOverlap.md",
            "api/Truncations.md"
        ]
    ]
)

# enables doc site deployment to Github Pages
deploydocs(
    repo="github.com/MSRudolph/PauliPropagation.jl.git",
)
