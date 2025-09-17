module ABM

using Agents
using Random
using Statistics

using ..Parameters: ModelParams, num_steps

const DEFAULT_POSITION = (0.0,)

@agent OpinionAgent ContinuousAgent{1} begin
    x::Float64
    r::Float64
    u::Float64
end

"""Create an `AgentBasedModel` initialised with opinion agents."""
function initialize_model(params::ModelParams)
    rng = Random.MersenneTwister(params.seed)
    space = ContinuousSpace(1; periodic=false)
    properties = Dict{Symbol, Any}(
        :params => params,
        :G => 0.0,
        :step => 0,
    )
    model = AgentBasedModel(OpinionAgent, space; properties=properties, rng=rng)

    for id in 1:params.N
        x = randn(rng)
        r = randn(rng)
        agent = OpinionAgent(id, DEFAULT_POSITION, x, r, x - r)
        add_agent!(agent, model)
    end

    model.G = isempty(model.agents) ? 0.0 : mean(agent.u for agent in allagents(model))
    return model
end

"""Single agent update according to the stochastic relaxation dynamics."""
function agent_step!(agent::OpinionAgent, model)
    params = model.params
    dt = params.dt
    rng = model.rng

    drift = (-params.λ * agent.x + params.κ * model.G) * dt
    diffusion = params.σ * sqrt(dt) * randn(rng)
    new_x = agent.x + drift + diffusion

    if abs(new_x - agent.r) >= params.Θ && rand(rng) < params.ϖ * dt
        new_x = agent.r
    end

    agent.x = new_x
    agent.u = agent.x - agent.r
    return
end

"""Update the cached mean-field order parameter stored on the model."""
function model_step!(model)
    model.step += 1
    model.G = isempty(model.agents) ? 0.0 : mean(agent.u for agent in allagents(model))
    return
end

"""Convert raw model data collected by `Agents.run!` into dictionaries."""
function _normalise_model_data(model_data, observables::Dict{Symbol, F}) where {F}
    if model_data === nothing
        return Dict{Symbol, Vector{Float64}}()
    end

    if hasmethod(Agents.collected_data, Tuple{typeof(model_data)})
        collected = Agents.collected_data(model_data)
        if collected !== nothing
            return _normalise_model_data(collected, observables)
        end
    end

    if model_data isa Dict
        return Dict(Symbol(k) => collect(v) for (k, v) in model_data)
    end

    if Base.hasproperty(model_data, :columnnames)
        columns = Dict{Symbol, Vector{Float64}}()
        for sym in keys(observables)
            if sym in propertynames(model_data)
                columns[sym] = collect(getproperty(model_data, sym))
            end
        end
        return columns
    end

    if model_data isa AbstractVector
        if !isempty(model_data) && first(model_data) isa NamedTuple
            rows = model_data
            return Dict{Symbol, Vector{Float64}}(
                sym => [row[sym] for row in rows if haskey(row, sym)]
                for sym in keys(observables)
            )
        elseif length(observables) == 1
            sym = first(keys(observables))
            return Dict(sym => collect(model_data))
        end
    end

    return Dict{Symbol, Vector{Float64}}()
end

function _default_observables()
    return Dict{Symbol, Function}(
        :G => model -> model.G,
        :mean_x => model -> (isempty(model.agents) ? 0.0 : mean(agent.x for agent in allagents(model))),
        :variance_u => model -> (isempty(model.agents) ? 0.0 : var(agent.u for agent in allagents(model))),
    )
end

function _coerce_observables(observables)
    if observables === nothing
        return _default_observables()
    elseif observables isa Dict
        return Dict(Symbol(k) => v for (k, v) in observables)
    else
        return Dict{Symbol, Function}(Symbol(name) => func for (name, func) in observables)
    end
end

"""
    run_abm(params; steps, observables, store)

Run the agent-based model for a number of steps while collecting observables.
Returns a named tuple with the evolved model, the raw result returned by
`Agents.run!`, and a dictionary of observable streams that is convenient for
analysis and plotting.
"""
function run_abm(params::ModelParams; steps::Union{Int,Nothing}=nothing,
                 observables=nothing, store=:model)
    nsteps = isnothing(steps) ? num_steps(params) : steps
    model_obs = _coerce_observables(observables)

    model = initialize_model(params)
    run_result = nothing
    try
        run_result = Agents.run!(model, agent_step!, model_step!, nsteps;
                                 mdata=model_obs, store=store)
    catch err
        run_result = err
    end

    if run_result isa Tuple
        agent_data, model_data = run_result
        data = _normalise_model_data(model_data, model_obs)
        return (model=model, agent_data=agent_data, model_data=model_data,
                data=data, observables=model_obs)
    elseif !(run_result isa Exception)
        data = _normalise_model_data(run_result, model_obs)
        return (model=model, agent_data=nothing, model_data=run_result,
                data=data, observables=model_obs)
    end

    model = initialize_model(params)
    collector = Agents.DataCollector(; model=model_obs, store=store)
    data = Dict{Symbol, Vector{Float64}}(sym => Float64[] for sym in keys(model_obs))

    for _ in 1:nsteps
        Agents.step!(model, agent_step!, model_step!, 1)
        try
            Agents.collect!(collector, model)
        catch
        end
        for (sym, obs) in model_obs
            push!(data[sym], obs(model))
        end
    end

    collected = try
        Agents.collected_data(collector)
    catch
        nothing
    end

    return (model=model, agent_data=nothing, model_data=collected,
            data=data, observables=model_obs)
end

export OpinionAgent, initialize_model, agent_step!, model_step!, run_abm

end # module ABM
