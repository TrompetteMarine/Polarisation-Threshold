module Analysis

using Statistics
using ..Core: Agent, ModelParams

function estimate_Vstar(agents::Vector{Agent})
    return var(agent.u for agent in agents)
end

function estimate_kappa_star(params::ModelParams, agents::Vector{Agent})
    V_star = estimate_Vstar(agents)
    return params.σ^2 / (2 * V_star)
end

function polarization_metric(x_history)
    return var(x_history[:, end])
end

export estimate_Vstar, estimate_kappa_star, polarization_metric

end # module Analysis
