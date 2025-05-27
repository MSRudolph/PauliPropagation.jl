# Generates HTML documentation from the contents of
# the docs folder. To generate, we must first setup
# a symlink from the repo README.md to src/index.md,
# in order re-use the README in Documenter.jl
# From the docs/ directory (containing this file):
#     cd src
#     ln -s ../../README.md index.md
#     cd ../
#
# This need only be done once per-machine. Then,
# generating/updating the doc is triggered via
#     julia --project make.jl
# 
# If triggered within a Github Action, the generated
# HTML files will then be committed to the 'gh-pages'
# branch, which Github Pages can be configured to
# display at msrudolph.github.io/PauliPropagation.jl/
# 
# Note documentation generated from non-main branches
# will be uploaded to subdomain /dev/, even when not
# from the 'dev' branch, and doc generated from pull
# requests will be uploaded to /previews/PR#.
#
# If this breaks, break Tyson's legs


using Documenter, PauliPropagation


# Generate doc HTML files, saved to build/
makedocs(
    # Add favicon.ico
    format=Documenter.HTML(
        assets=[
            "assets/favicon.ico",
        ],
    ), sitename="PauliPropagation.jl",

    # determines site layout
    pages=[

        # index.md does not exist; it is a symlink
        # to the repo's README.md file, created as
        # per the comments above, to avoid duplicating
        # the README.md contents into Documenter.jl 
        # pages. We manually override its name in the
        # left navbar to be "Introduction"
        "Home" => "index.md",

        # these other 'top-level' files DO exist, and
        # have names inferred from their section names

        # TODO: add this back once we know how to embed the Jupyter notebooks
        # "tutorials.md",

        # these 'lower-level' files also exist, and will
        # be grouped under an 'API' section in the navbar
        "API" => [
            "api/PauliAlgebra.md",
            "api/Gates.md",
            "api/Circuits.md",
            "api/Propagation.md",
            "api/StateOverlap.md",
            "api/PathProperties.md",
            "api/PauliTransferMatrix.md",
            "api/Surrogate.md",
            "api/NumericalCertificates.md",
            "api/Truncations.md"
        ]
    ]
)


# When run from a Github Action, commit those files to the 'gh-pages' branch.
# If the action's invoking branch is not main, the files are uploaded under /dev/
deploydocs(
    repo="github.com/MSRudolph/PauliPropagation.jl.git",

    # Enable generation of doc from PRs, under a /previews/PR## sub-domain.
    # Beware that this requires the Github Action was explicitly triggered by
    # a 'pull_request' event (not a 'push')
    push_preview=true,

    # Specify the name of our develop branch (Documenter.jl seems unable to 
    # auto-infer it) so that changes thereto generate preview-documentation 
    devbranch="dev"
)


# Once 'gh-pages' branch is updated, and Github Pages has been configured to
# publish files from that branch, the documentation is visible at either:
# - msrudolph.github.io/PauliPropagation.jl/
# - msrudolph.github.io/PauliPropagation.jl/dev/
# - msrudolph.github.io/PauliPropagation.jl/previews/PR#
# where # above is replaced with the pull request number.
#
# These "doc clones" are deleted whenever a commit is pushed to the main
# branch (signifying a version release), so that development history does
# not bloat the repo. Deletion is performed by the 'tidy-doc' CI job.
