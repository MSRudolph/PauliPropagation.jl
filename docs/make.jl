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
    devbranch = "doc-ci",
    repo = "github.com/MSRudolph/PauliPropagation.jl.git",
)


# TODO!
# change "doc-ci" above to "main" before merging PR to main
