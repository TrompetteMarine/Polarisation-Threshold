using Statistics
using DifferentialEquations
using ..Core: ModelParams, Agent, initialize_agents, ou_drift!, ou_diffusion!
import ..Simulation: reset_condition, reset_affect!

function compute_mean_field(agents::Vector{Agent})
    return mean(agent.x - agent.r for agent in agents)
end

function simulate!(agents::Vector{Agent}, params::ModelParams)
    times = 0:params.dt:params.T
    G_history = zeros(length(times))
    x_history = zeros(params.N, length(times))
    reset_cb = DiscreteCallback(reset_condition, reset_affect!)

    for (t_idx, t) in enumerate(times)
        G = compute_mean_field(agents)
        G_history[t_idx] = G
        for (i, agent) in enumerate(agents)
            p = (λ=params.λ, κ=params.κ, G=G, σ=params.σ, Θ=params.Θ, ϖ=params.ϖ)
            prob = SDEProblem(ou_drift!, ou_diffusion!, [agent.x, agent.r], (t, t + params.dt), p)
            sol = solve(prob, EM(), dt=params.dt, callback=reset_cb)
            agent.x = sol.u[end][1]
            agent.u = agent.x - agent.r
            x_history[i, t_idx] = agent.x
        end
    end
    return (times=times, G_history=G_history, x_history=x_history)
end
