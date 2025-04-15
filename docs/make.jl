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

    # do not place in-development doc under a /dev/ sub-domain,
    # since we'll instead preview doc straight from PRs (as below)
    # and do not wish to encourage distributing the wrong URL
    devurl="",

    # enable generation of doc from PRs, under a /previews/PR## sub-domain
    push_preview=true,

    # DEBUG
    # this forces deploydocs to upload to gh-pages when the invoking branch is
    # not "main" (the default). Otherwise, despite push_preview=true, the upload
    # is cancelled because ENV["GITHUB_REF"] does not match devbranch="main". Bug??
    devbranch="doc-fix"
)
