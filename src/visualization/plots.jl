module Visualization

using Statistics
using Plots
using ..ModelCore: ModelParams
using ..Simulation: simulate!, initialize_agents

function bifurcation_plot(params::ModelParams; κ_range=0.0:0.01:2.0)
    bifurcation_data = Float64[]
    for κ in κ_range
        step_params = ModelParams(
            λ = params.λ,
            σ = params.σ,
            κ = κ,
            Θ = params.Θ,
            ϖ = params.ϖ,
            N = params.N,
            T = params.T,
            dt = params.dt,
            seed = params.seed,
        )
        agents = initialize_agents(step_params)
        times, G_history, _ = simulate!(agents, step_params)
        eq_G = mean(abs.(G_history[end÷10:end]))
        push!(bifurcation_data, eq_G)
    end
    plot(κ_range, bifurcation_data,
        xlabel="Peer Coupling Strength (κ)",
        ylabel="Equilibrium Polarization (|G|)",
        title="Bifurcation Diagram",
        label="Polarization",
        linewidth=2)
end

function phase_portrait(agents, params::ModelParams)
    scatter([agent.r for agent in agents], [agent.x for agent in agents],
            xlabel="Anchor (r)",
            ylabel="Belief (x)",
            title="Phase Portrait (κ = $(params.κ))",
            label="Agents",
            alpha=0.5)
    plot!([-3, 3], [-3, 3], linestyle=:dash, label="x = r")
end

function plot_mean_field(times, G_history, params::ModelParams)
    plot(times, G_history,
         xlabel="Time",
         ylabel="Mean-Field Term (G)",
         title="Evolution of Mean-Field Term (κ = $(params.κ))",
         label="G(t)",
         linewidth=2)
end

export bifurcation_plot, phase_portrait, plot_mean_field

end # module Visualization
