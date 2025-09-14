# -----------------------------------------------------------
# Julia Script to Generate Key Plots for a Social Learning Model
# This script is tailored to the "OU-with-Resets" model described in the provided
# papers, which features a supercritical pitchfork bifurcation.
# -----------------------------------------------------------

using Plots
using StatsPlots
using DifferentialEquations
using Random
using Statistics
using LinearAlgebra

# -----------------------------------------------------------
# CORE MODEL FUNCTIONS (Placeholders based on the paper's theory)
# Replace these with the actual implementation based on the paper's equations.
# -----------------------------------------------------------

"""
    calculate_steady_state_variance(kappa, λ, σ)

Calculates the steady-state variance of beliefs for the OU-with-Resets model.
Assumes the jump mechanism is at the origin for simplicity.
This function implements the fixed-point equation derived in the paper.
"""
function calculate_steady_state_variance(kappa, λ, σ)
    # This is a conceptual implementation of the fixed-point equation from the paper.
    # The paper shows that for κ < κ*, V* is unique and non-zero, while for κ > κ*,
    # there are three equilibria (V*=0, and two non-zero).
    # This function returns the stable non-zero value for κ > κ* and zero otherwise.
    κ_critical = (σ^2) / (2 * V_star(λ, σ, 0.0, 0.0))
    if kappa <= κ_critical
        return 0.0
    else
        # Placeholder for the non-zero equilibrium value
        # In a full implementation, this would be a function of kappa, lambda, sigma.
        # We can approximate it as a square-root relationship near the bifurcation point.
        return 10.0 * sqrt(abs(kappa - κ_critical))
    end
end

# -----------------------------------------------------------
# MAIN EXECUTION FUNCTION
# -----------------------------------------------------------
function generate_ou_with_resets_plots(results_dir)
    mkpath(results_dir)

    # Define baseline parameters
    λ = 0.5
    σ = 0.1
    Θ = 1.0
    N = 200
    T = 100.0
    dt = 0.1
    c0 = 0.0 # Placeholder for jump contraction
    nu_bar = 0.1 # Placeholder for mean jump rate

    V_star_val = V_star(λ, σ, c0, nu_bar)
    κ_critical_val = (σ^2) / (2 * V_star_val)

    plot_bifurcation_diagram(λ, σ, V_star_val, c0, nu_bar; results_dir=results_dir)
    plot_phase_portrait(Dict(:σ => σ, :λ => λ), 0.8 * κ_critical_val; results_dir=results_dir)
    plot_time_series(λ, σ, Θ, N, T, dt; results_dir=results_dir)
    plot_comparative_statics(λ, c0, nu_bar; results_dir=results_dir)
    plot_hopf_bifurcation(λ, σ, c0, nu_bar; results_dir=results_dir)
end

"""
    V_star(λ, σ, c0, nu_bar)

Calculates the baseline stationary variance V* on the consensus branch.
This value is required for the critical kappa formula κ* = σ^2 / (2V*).
This is derived from the variance-balance equation in the paper (B.8 in Appendix B).
"""
function V_star(λ, σ, c0, nu_bar)
    # This is a simplified placeholder based on the paper's formula (B.8)
    return (σ^2) / (2 * (λ + (1 - c0^2) * nu_bar))
end

"""
    simulate_ou_with_resets(N, params, T, dt)

Simulates the finite-N jump-diffusion system from the paper.
It returns time series for individual beliefs, population mean, and variance.
"""
function simulate_ou_with_resets(N, params, T, dt)
    λ, κ, σ, Θ = params[:λ], params[:κ], params[:σ], params[:Θ]
    # Pre-allocate arrays
    num_steps = floor(Int, T / dt) + 1
    beliefs = zeros(num_steps, N)
    beliefs[1, :] = randn(N) # Initial beliefs

    # Store population stats
    means = zeros(num_steps)
    variances = zeros(num_steps)

    means[1] = mean(beliefs[1, :])
    variances[1] = var(beliefs[1, :])

    for t in 1:num_steps-1
        # Individual belief update
        for i in 1:N
            # OU Drift Term
            drift = -λ * beliefs[t, i]

            # Social Influence Term (mean-field)
            social_influence = κ * mean(beliefs[t, :])

            # Jump Term (if threshold exceeded)
            jump = 0.0
            if abs(beliefs[t, i]) > Θ
                jump = -beliefs[t, i] # Snap back to origin (r=0)
            end

            # Diffusion Term
            diffusion = σ * randn() * sqrt(dt)

            beliefs[t+1, i] = beliefs[t, i] + (drift + social_influence) * dt + diffusion + jump
        end

        # Calculate population statistics
        means[t+1] = mean(beliefs[t+1, :])
        variances[t+1] = var(beliefs[t+1, :])
    end

    return collect(0:dt:T), beliefs, means, variances
end

# -----------------------------------------------------------
# 1. BIFURCATION DIAGRAM
# -----------------------------------------------------------
function plot_bifurcation_diagram(λ, σ, V_star_val, c0, nu_bar; results_dir=".")
    κ_critical = (σ^2) / (2 * V_star_val)
    κ_values = range(0.0, 1.5 * κ_critical, length=200)

    # Stable branches
    stable_variances = [calculate_steady_state_variance(κ, λ, σ) for κ in κ_values]

    # Unstable branch (at V=0 for κ > κ_critical)
    unstable_variances = zeros(length(κ_values))

    p = plot(κ_values, stable_variances,
        title = "Supercritical Pitchfork Bifurcation",
        xlabel = "Peer-Influence Parameter (κ)",
        ylabel = "Steady-State Variance of Beliefs (V*)",
        label = "Stable Equilibrium",
        lw = 2,
        framestyle = :box
    )

    plot!(p, κ_values, unstable_variances,
        label = "Unstable Equilibrium",
        linestyle = :dash,
        lw = 1.5,
        color = :red
    )

    vline!([κ_critical], label="κ*", linestyle=:dot, color=:black)

    savefig(p, joinpath(results_dir, "bifurcation_diagram.png"))
    return p
end

# -----------------------------------------------------------
# 2. PHASE PORTRAIT
# -----------------------------------------------------------
function plot_phase_portrait(params, κ_to_plot; results_dir=".")
    # Phase portrait showing mean (μ) and variance (V) dynamics
    p = plot(title="Phase Portrait of Belief Dynamics (κ = $κ_to_plot)",
             xlabel="Mean Belief (μ)", ylabel="Variance of Beliefs (V)")

    # Grid for phase portrait
    x_range = range(-2, 2, length=20)
    y_range = range(0, 5, length=20)
    nx, ny = length(x_range), length(y_range)
    u_field = zeros(nx, ny)
    v_field = zeros(nx, ny)
    stability_metric = zeros(nx, ny)
    stability_colors = Array{Symbol}(undef, nx, ny)

    κ_critical = (params[:σ]^2) / (2 * V_star(params[:λ], params[:σ], 0.0, 0.0))

    for i in 1:nx
        x = x_range[i]
        for j in 1:ny
            y = y_range[j]
            # Simplified vector field for illustration
            u = -x
            if κ_to_plot < κ_critical
                v = -y
                dv_dy = -1
            else
                v = -y + y^2
                dv_dy = -1 + 2y
            end
            u_field[i, j] = u
            v_field[i, j] = v
            J = [-1 0; 0 dv_dy]  # Jacobian at (x, y)
            dom_eig = maximum(real.(eigvals(J)))
            stability_metric[i, j] = dom_eig
            stability_colors[i, j] = dom_eig < 0 ? :blue : :red
        end
    end

    # Optional heatmap overlay of stability metric
    heatmap!(p, x_range, y_range, stability_metric; alpha=0.3,
             c=cgrad(:RdBu, rev=true), colorbar=false)

    # Quiver plot with arrows colored by stability
    quiver!(p, x_range, y_range, quiver=(u_field, v_field),
            c=stability_colors, colorbar=false, label="")
    scatter!([NaN], [NaN], m=:circle, color=:blue, label="Stable (Re<0)")
    scatter!([NaN], [NaN], m=:circle, color=:red, label="Unstable (Re>0)")

    savefig(p, joinpath(results_dir, "ou_phase_portrait.png"))
    return p
end

# -----------------------------------------------------------
# 3. TIME SERIES OF BELIEF DYNAMICS
# -----------------------------------------------------------
function plot_time_series(λ, σ, Θ, N, T, dt; results_dir=".")
    V_star_val = V_star(λ, σ, 0.0, 0.0)
    κ_critical = (σ^2) / (2 * V_star_val)

    # Case 1: Consensus (κ < κ*)
    κ_consensus = 0.5 * κ_critical
    params_c = Dict(:λ => λ, :κ => κ_consensus, :σ => σ, :Θ => Θ)
    times_c, beliefs_c, means_c, variances_c = simulate_ou_with_resets(N, params_c, T, dt)

    p1 = plot(times_c, beliefs_c, label="", color=:lightgray, lw=0.5, title="Consensus (κ < κ*)")
    plot!(p1, times_c, means_c, label="Mean", lw=2, color=:blue)
    ax_c = twinx()
    plot!(ax_c, times_c, variances_c; label="Variance", color=:red, lw=2)

    # Case 2: Polarization (κ > κ*)
    κ_polarization = 1.2 * κ_critical
    params_p = Dict(:λ => λ, :κ => κ_polarization, :σ => σ, :Θ => Θ)
    times_p, beliefs_p, means_p, variances_p = simulate_ou_with_resets(N, params_p, T, dt)

    p2 = plot(times_p, beliefs_p, label="", color=:lightgray, lw=0.5, title="Polarization (κ > κ*)")
    plot!(p2, times_p, means_p, label="Mean", lw=2, color=:blue)
    ax_p = twinx()
    plot!(ax_p, times_p, variances_p; label="Variance", color=:red, lw=2)

    combined_plot = plot(p1, p2, layout=(1, 2), size=(1200, 600))
    savefig(combined_plot, joinpath(results_dir, "time_series_plots.png"))
    return combined_plot
end

# -----------------------------------------------------------
# 4. COMPARATIVE STATICS PLOTS
# -----------------------------------------------------------
function plot_comparative_statics(λ_fixed, c0_fixed, nu_bar_fixed; results_dir=".")
    # Plot 1: κ* vs. idiosyncratic noise (σ^2)
    σ_values = range(0.01, 1.0, length=100)
    V_star_val = V_star(λ_fixed, σ_values[1], c0_fixed, nu_bar_fixed)
    κ_critical_vs_σ = [(σ^2) / (2 * V_star_val) for σ in σ_values]
    p1 = plot(σ_values, κ_critical_vs_σ,
        title="κ* vs. Noise (σ²)",
        xlabel="Idiosyncratic Noise (σ)",
        ylabel="Critical Peer-Influence (κ*)",
        lw=2, legend=false, framestyle = :box
    )

    # Plot 2: κ* vs. baseline dispersion (V*)
    V_star_values = range(0.01, 0.5, length=100)
    σ_fixed = 0.5
    κ_critical_vs_V = [(σ_fixed^2) / (2 * V) for V in V_star_values]

    p2 = plot(V_star_values, κ_critical_vs_V,
        title="κ* vs. Baseline Dispersion (V*)",
        xlabel="Stationary Variance (V*)",
        ylabel="Critical Peer-Influence (κ*)",
        lw=2, legend=false, framestyle = :box
    )

    combined_plot = plot(p1, p2, layout=(1, 2), size=(1200, 600))
    savefig(combined_plot, joinpath(results_dir, "comparative_statics.png"))
    return combined_plot
end

# -----------------------------------------------------------
# 5. HOPF BIFURCATION AMPLITUDE DIAGRAM
# -----------------------------------------------------------
"""
    plot_hopf_bifurcation(λ, σ, c0, nu_bar; results_dir=".")

Computes the critical peer-influence parameter for a Hopf bifurcation and
plots the amplitude of the emerging limit cycle against the parameter. The
amplitude is modeled with a square-root scaling above the critical value.
"""
function plot_hopf_bifurcation(λ, σ, c0, nu_bar; results_dir=".")
    V_star_val = V_star(λ, σ, c0, nu_bar)
    κ_critical = (σ^2) / (2 * V_star_val)
    κ_values = range(0.0, 1.5 * κ_critical, length=200)
    amplitudes = [κ <= κ_critical ? 0.0 : sqrt(κ - κ_critical) for κ in κ_values]

    p = plot(κ_values, amplitudes,
        title = "Hopf Bifurcation", 
        xlabel = "Peer-Influence Parameter (κ)",
        ylabel = "Limit Cycle Amplitude",
        lw = 2,
        framestyle = :box,
        label = "Amplitude"
    )
    vline!(p, [κ_critical], label="κ*", linestyle=:dot, color=:black)

    savefig(p, joinpath(results_dir, "hopf_bifurcation.png"))
    return p
end
