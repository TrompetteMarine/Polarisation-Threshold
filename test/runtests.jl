using Test
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using PolarisationThreshold
using DynamicalSystems: ContinuousDynamicalSystem

include("abm_tests.jl")
include("meanfield_tests.jl")
include("analysis_tests.jl")
