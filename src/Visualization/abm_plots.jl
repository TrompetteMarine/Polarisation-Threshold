const DEFAULT_RESULTS_DIR = joinpath(@__DIR__, "..", "..", "results")

function _ensure_results_dir(results_dir)
    path = abspath(results_dir === nothing ? DEFAULT_RESULTS_DIR : results_dir)
    mkpath(path)
    return path
end

"""Plot a collection of observables from the ABM data dictionary."""
function plot_model_observables(model_data; selection=nothing, results_dir=DEFAULT_RESULTS_DIR, filename="abm_observables.png")
    observables = selection === nothing ? collect(Base.keys(model_data)) : selection
    plt = plot()
    for key in observables
        sym = Symbol(key)
        if haskey(model_data, sym)
            values = collect(model_data[sym])
            plot!(plt, values; label=String(sym))
        end
    end
    xlabel!(plt, "Step")
    ylabel!(plt, "Value")
    title!(plt, "ABM Observables")
    path = joinpath(_ensure_results_dir(results_dir), filename)
    savefig(plt, path)
    return plt
end

"""Plot the mean-field series from ABM output."""
function plot_mean_field_series(model_data; results_dir=DEFAULT_RESULTS_DIR, filename="abm_mean_field.png")
    if !haskey(model_data, :G)
        error("Observable :G not found in model data")
    end
    values = collect(model_data[:G])
    plt = plot(values; xlabel="Step", ylabel="Mean field G", label="G", title="Mean field trajectory")
    path = joinpath(_ensure_results_dir(results_dir), filename)
    savefig(plt, path)
    return plt
end

"""Plot the variance series captured during ABM simulation."""
function plot_variance_series(model_data; results_dir=DEFAULT_RESULTS_DIR, filename="abm_variance.png")
    sym = :variance_u
    if !haskey(model_data, sym)
        error("Observable :variance_u not found in model data")
    end
    values = collect(model_data[sym])
    plt = plot(values; xlabel="Step", ylabel="Variance", label="var(u)", title="Opinion variance")
    path = joinpath(_ensure_results_dir(results_dir), filename)
    savefig(plt, path)
    return plt
end

