"""Plot the time evolution of the first component of an orbit."""
function plot_orbit_trajectory(orbit; results_dir=DEFAULT_RESULTS_DIR, filename="mean_field_orbit.png")
    times = get(orbit, :times, nothing)
    states = get(orbit, :states, nothing)
    if times === nothing || states === nothing
        error("Orbit data must contain :times and :states entries")
    end
    values = states[1, :]
    plt = plot(times, values; xlabel="Time", ylabel="G", label="G(t)", title="Mean-field orbit")
    path = joinpath(_ensure_results_dir(results_dir), filename)
    savefig(plt, path)
    return plt
end

"""Render a bifurcation diagram from aggregated data."""
function plot_bifurcation_diagram(data; results_dir=DEFAULT_RESULTS_DIR, filename="bifurcation.png", alpha=0.4)
    params = get(data, :parameters, nothing)
    observable = get(data, :observable, nothing)
    if params === nothing || observable === nothing
        error("Bifurcation data must contain :parameters and :observable entries")
    end
    plt = scatter(params, observable; xlabel="Parameter", ylabel="Observable", alpha=alpha, markersize=2, label="orbit")
    path = joinpath(_ensure_results_dir(results_dir), filename)
    savefig(plt, path)
    return plt
end
