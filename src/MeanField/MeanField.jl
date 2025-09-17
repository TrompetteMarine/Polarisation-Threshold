module MeanField

using DynamicalSystems

using ..Parameters: ModelParams

"""Construct the continuous mean-field dynamical system for a given parameter set."""
function mean_field_system(params::ModelParams; initial_state=nothing, reduced::Bool=false)
    function flow!(du, u, p, t)
        x = u[1]
        r = u[2]
        G = x - r
        drift = -params.λ * x + params.κ * G
        if reduced && abs(G) >= params.Θ
            drift -= params.ϖ * (G - sign(G) * params.Θ)
        end
        du[1] = drift
        du[2] = 0.0
    end

    u0 = initial_state === nothing ? [0.0, 0.0] : Float64.(collect(initial_state))
    return ContinuousDynamicalSystem(flow!, u0, params)
end

"""Construct a reduced one-dimensional flow for the mean-field order parameter."""
function reduced_reset_map(params::ModelParams; initial_state=nothing)
    function flow!(du, u, p, t)
        G = u[1]
        drift = (-params.λ + params.κ) * G
        if abs(G) >= params.Θ
            drift -= params.ϖ * (G - sign(G) * params.Θ)
        end
        du[1] = drift
    end

    u0 = initial_state === nothing ? [0.0] : Float64.(collect(initial_state))
    return ContinuousDynamicalSystem(flow!, u0, params)
end

function _default_dt(params)
    hasproperty(params, :dt) ? params.dt : 1.0
end

"""Generate orbit data from a dynamical system."""
function orbit_data(system::ContinuousDynamicalSystem, initial_state=nothing; steps::Int=1000, discard::Int=0, dt=nothing)
    if initial_state !== nothing
        set_state!(system, Float64.(collect(initial_state)))
    end
    Δt = dt === nothing ? _default_dt(system.p) : dt
    traj = trajectory(system, steps + discard; Δt=Δt)
    data = Array(traj)
    if discard > 0
        data = data[:, discard+1:end]
    end
    times = collect(0:Δt:Δt * (size(data, 2) - 1))
    return (times=times, states=data)
end

"""Sweep a parameter range and collect observable values for bifurcation analysis."""
function bifurcation_data(system_factory, param_range; observable = state -> state[1], steps::Int=1000, discard::Int=200, dt=nothing)
    param_values = Float64[]
    observable_values = Float64[]
    states_by_param = Dict{Float64, Matrix{Float64}}()

    for pval in param_range
        system = system_factory(pval)
        data = orbit_data(system; steps=steps, discard=discard, dt=dt)
        states = data.states
        states_by_param[Float64(pval)] = states
        for col in axes(states, 2)
            push!(param_values, Float64(pval))
            push!(observable_values, observable(view(states, :, col)))
        end
    end

    return (parameters=param_values, observable=observable_values, states=states_by_param)
end

export mean_field_system, reduced_reset_map, orbit_data, bifurcation_data

end # module MeanField
