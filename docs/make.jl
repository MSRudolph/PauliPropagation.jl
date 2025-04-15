using Documenter, PauliPropagation

# determines doc site layout
makedocs(
    sitename="Pauli Propagation",
    pages=[
        "introduction.md",
        "installation.md",
        "tutorials.md",
        "index.md",
    ]
)

# enables doc site deployment to Github Pages
deploydocs(
    repo="github.com/MSRudolph/PauliPropagation.jl.git",
)
