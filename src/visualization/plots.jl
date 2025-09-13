module Visualization

using Statistics
using Plots
using ..ModelCore: ModelParams
using ..Simulation: simulate!, initialize_agents

"""
    moving_average(data, window)

Simple moving average used to smooth the mean orbit line in the bifurcation
diagram. Points containing `NaN` are skipped when computing the local mean.
"""
function moving_average(data::AbstractVector{<:Real}, window::Int)
    n = length(data)
    smoothed = Vector{Float64}(undef, n)
    for i in 1:n
        lo = max(1, i - window ÷ 2)
        hi = min(n, i + window ÷ 2)
        vals = data[lo:hi]
        vals = filter(!isnan, vals)
        smoothed[i] = isempty(vals) ? NaN : mean(vals)
    end
    return smoothed
end

"""
    bifurcation_plot(params; κ_range)

Generate a bifurcation diagram by sweeping over `κ_range`. For each `κ`, the
system is simulated, an initial transient (first 20% of steps) is discarded and
the remaining orbit values of the mean-field term `G` are collected. Points are
classified as stable when `abs(G) ≤ params.Θ` and unstable otherwise. Stable and
unstable orbits are plotted in different colours with low opacity. A smoothed
line through the mean of the stable orbits is overlaid for reference.
"""
function bifurcation_plot(params::ModelParams; κ_range=0.0:0.01:2.0)
    stable_κ = Float64[]
    stable_G = Float64[]
    unstable_κ = Float64[]
    unstable_G = Float64[]
    stable_map = Dict{Float64, Vector{Float64}}()

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

        transient_idx = floor(Int, 0.2 * length(G_history)) + 1
        G_vals = G_history[transient_idx:end]

        for G in G_vals
            if abs(G) <= params.Θ
                push!(stable_κ, κ)
                push!(stable_G, G)
                push!(get!(stable_map, κ, Float64[]), G)
            else
                push!(unstable_κ, κ)
                push!(unstable_G, G)
            end
        end
    end

    plt = scatter(stable_κ, stable_G,
                  xlabel="Peer Coupling Strength (κ)",
                  ylabel="Mean-Field Term (G)",
                  title="Bifurcation Diagram",
                  color=:blue, alpha=0.3, markersize=2, label="Stable")
    scatter!(plt, unstable_κ, unstable_G,
             color=:red, alpha=0.3, markersize=2, label="Unstable")

    stable_means = [haskey(stable_map, κ) ? mean(stable_map[κ]) : NaN for κ in κ_range]
    smoothed = moving_average(stable_means, 5)
    valid = .!isnan(smoothed)
    plot!(plt, collect(κ_range)[valid], smoothed[valid],
          color=:blue, linewidth=2, label="Stable mean")

    return plt
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
