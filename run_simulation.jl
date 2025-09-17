#!/usr/bin/env julia
using Pkg
Pkg.activate(@__DIR__)
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

using PolarisationThreshold

function main()
    params = ModelParams()

    abm_result = ABM.run_abm(params)
    Visualization.plot_mean_field_series(abm_result.data)
    Visualization.plot_variance_series(abm_result.data)

    system = MeanField.mean_field_system(params)
    orbit = MeanField.orbit_data(system; steps=1500, discard=200)
    Visualization.plot_orbit_trajectory(orbit)

    factory = κ -> MeanField.mean_field_system(update(params; κ=κ); reduced=true)
    bif_data = MeanField.bifurcation_data(factory, 0.0:0.05:2.0; observable = state -> state[1], steps=1200, discard=200)
    Visualization.plot_bifurcation_diagram(bif_data)

    thresholds = Analysis.compare_thresholds(abm_result.data, orbit, params)
    println("κ* estimates: ", thresholds)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
