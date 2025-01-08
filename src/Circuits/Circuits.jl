### Circuits/Circuits.jl
##
# This module contains functions for building circuits for unitary evolution.
# Circuits are generally vectors of `Gate` objects.
# It provides functions to generate topologies, build circuits, and utility functions.
##
###

include("utils.jl")
include("topologies.jl")
include("builders.jl")