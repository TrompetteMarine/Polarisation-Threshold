using .Core: ModelParams, Agent, initialize_agents
using .Simulation: simulate!, compute_mean_field, reset_condition, reset_affect!
using .Visualization: bifurcation_plot, phase_portrait, plot_mean_field
using .Analysis: estimate_kappa_star, polarization_metric

function main()
    params = ModelParams(κ=0.5, N=100, T=50.0)
    agents = initialize_agents(params)
    times, G_history, x_history = simulate!(agents, params)

    bifurcation_plot(params, κ_range=0.0:0.01:1.0)
    phase_portrait(agents, params)
    plot_mean_field(times, G_history, params)

    κ_star = estimate_kappa_star(params, agents)
    println("Estimated critical coupling κ*: ", κ_star)
    println("Polarization metric: ", polarization_metric(x_history))
end

main()
