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
    G = mean(agent.x - agent.r for agent in agents)

    r_vals = range(-3, 3, length=20)
    x_vals = range(-3, 3, length=20)
    r_grid = repeat(collect(r_vals)', length(x_vals), 1)
    x_grid = repeat(collect(x_vals), 1, length(r_vals))

    dx = -params.λ .* x_grid .+ params.κ * G
    dr = zeros(size(dx))

    stable = abs.(x_grid .- r_grid) .<= params.Θ
    stable_int = Int.(stable)

    heatmap(r_vals, x_vals, stable_int;
            c=[:red, :green],
            alpha=0.3,
            xlabel="Anchor (r)",
            ylabel="Belief (x)",
            title="Phase Portrait (κ = $(params.κ))",
            legend=false)
    quiver!(r_grid, x_grid, dr, dx; color=:black, legend=false)

    agent_r = [agent.r for agent in agents]
    agent_x = [agent.x for agent in agents]
    agent_stable = abs.(agent_x .- agent_r) .<= params.Θ
    agent_colors = [s ? :green : :red for s in agent_stable]

    scatter!(agent_r, agent_x;
             color=agent_colors,
             label="Agents",
             alpha=0.8)
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
