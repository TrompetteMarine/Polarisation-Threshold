#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using JLD2
using PolarisationThreshold

function main(; κ_range=0.0:0.05:2.0)
    params = ModelParams()
    system = MeanField.mean_field_system(params)
    orbit = MeanField.orbit_data(system; steps=2000, discard=200)

    results_dir = joinpath(@__DIR__, "..", "results")
    mkpath(results_dir)

    Visualization.plot_orbit_trajectory(orbit; results_dir=results_dir)

    factory = κ -> MeanField.mean_field_system(update(params; κ=κ); reduced=true)
    bif_data = MeanField.bifurcation_data(factory, κ_range; observable = state -> state[1], steps=1200, discard=200)
    Visualization.plot_bifurcation_diagram(bif_data; results_dir=results_dir)

    @save joinpath(results_dir, "mean_field_results.jld2") params orbit bifurcation=bif_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
