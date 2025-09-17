module Parameters

using Parameters: @with_kw

@with_kw struct ModelParams
    λ::Float64 = 1.0
    σ::Float64 = 0.5
    κ::Float64 = 0.15
    Θ::Float64 = 1.0
    ϖ::Float64 = 0.2
    N::Int = 100
    T::Float64 = 100.0
    dt::Float64 = 0.01
    seed::Int = 1234
end

"""Return the number of discrete simulation steps implied by `params`."""
function num_steps(params::ModelParams)
    return max(1, round(Int, params.T / params.dt))
end

"""Convert the parameter struct to a named tuple."""
function as_namedtuple(params::ModelParams)
    return (
        λ = params.λ,
        σ = params.σ,
        κ = params.κ,
        Θ = params.Θ,
        ϖ = params.ϖ,
        N = params.N,
        T = params.T,
        dt = params.dt,
        seed = params.seed,
    )
end

"""Create a modified copy of `params` with updated keyword arguments."""
function update(params::ModelParams; kwargs...)
    base = as_namedtuple(params)
    merged = merge(base, NamedTuple(kwargs))
    return ModelParams(; merged...)
end

export ModelParams, num_steps, as_namedtuple, update

end # module Parameters
