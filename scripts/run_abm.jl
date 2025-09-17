#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using JLD2
using PolarisationThreshold

function main(; steps=nothing)
    params = ModelParams()
    result = ABM.run_abm(params; steps=steps)
    results_dir = joinpath(@__DIR__, "..", "results")
    mkpath(results_dir)

    Visualization.plot_mean_field_series(result.data; results_dir=results_dir)
    Visualization.plot_variance_series(result.data; results_dir=results_dir)
    Visualization.plot_model_observables(result.data; results_dir=results_dir)

    @save joinpath(results_dir, "abm_results.jld2") params result_data=result.data observables=result.observables
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
