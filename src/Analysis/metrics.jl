using Statistics

using ..Parameters: ModelParams

"""Extract the variance stream from ABM model data."""
function variance_stream(model_data)
    for key in keys(model_data)
        if occursin("variance", String(key))
            return collect(model_data[key])
        end
    end
    error("No variance observable found in collected data")
end

"""Estimate the critical coupling κ* from ABM variance data."""
function kappa_star_from_abm(model_data, params::ModelParams; window::Int=10)
    variances = variance_stream(model_data)
    if isempty(variances)
        return NaN
    end
    tail = variances[max(1, end - window + 1):end]
    mean_variance = mean(tail)
    if mean_variance == 0
        return Inf
    end
    return params.σ^2 / (2 * mean_variance)
end

"""Compute a deterministic stability measure from mean-field trajectories."""
function stability_fraction(orbit, params::ModelParams)
    states = orbit.states
    if isempty(states)
        return 0.0
    end
    within = 0
    total = size(states, 2)
    for col in axes(states, 2)
        G = states[1, col]
        within += abs(G) <= params.Θ ? 1 : 0
    end
    return within / total
end

"""Estimate κ* using the variance of a deterministic orbit."""
function kappa_star_from_orbit(orbit, params::ModelParams)
    states = orbit.states
    if isempty(states)
        return NaN
    end
    values = vec(states[1, :])
    variance = var(values)
    if variance == 0
        return Inf
    end
    return params.σ^2 / (2 * variance)
end

"""Compare stochastic and deterministic κ* estimates."""
function compare_thresholds(abm_data, orbit, params::ModelParams; window::Int=10)
    κ_abm = kappa_star_from_abm(abm_data, params; window=window)
    κ_det = kappa_star_from_orbit(orbit, params)
    return (abm=κ_abm, deterministic=κ_det, difference=κ_abm - κ_det)
end

export variance_stream, kappa_star_from_abm, stability_fraction, kappa_star_from_orbit, compare_thresholds
