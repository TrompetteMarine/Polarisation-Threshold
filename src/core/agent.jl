using Random

mutable struct Agent
    x::Float64
    r::Float64
    u::Float64
end

function initialize_agents(params::ModelParams)
    Random.seed!(params.seed)
    agents = [Agent(randn(), randn(), 0.0) for _ in 1:params.N]
    for agent in agents
        agent.u = agent.x - agent.r
    end
    return agents
end
