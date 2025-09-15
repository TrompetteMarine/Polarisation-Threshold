# ===========================================================
# Plot generation for the OU-with-Resets toy model
# (self-contained; no external project code required)
#
# Run:
#   julia --project=. scripts/generate_plots.jl
# or just:
#   julia scripts/generate_plots.jl
# ===========================================================

using Plots
using Random
using Statistics
using LinearAlgebra

# -----------------------------------------------------------
# Utilities
# -----------------------------------------------------------

"""
    V_star(lambda, sigma, c0, nubar)

Toy stationary variance on the consensus branch (placeholder).
Used to define κ* = σ^2 / (2 V*). Keep c0, nubar for compatibility.
"""
function V_star(lambda::Float64, sigma::Float64, c0::Float64, nubar::Float64)
    # Simple OU-like balance with additional dissipation from jump mechanism
    return (sigma^2) / (2 * (lambda + (1 - c0^2) * nubar + eps()))  # eps() to avoid div-by-zero
end

"""
    calculate_steady_state_variance(kappa, lambda, sigma, c0, nubar; scale=10.0)

Toy bifurcation branch: zero below κ*, positive above with sqrt-scaling.
This is a stand-in visualization (not a calibrated formula).
"""
function calculate_steady_state_variance(kappa::Float64, lambda::Float64, sigma::Float64,
                                         c0::Float64, nubar::Float64; scale::Float64=10.0)
    kappa_crit = (sigma^2) / (2 * V_star(lambda, sigma, c0, nubar))
    if kappa <= kappa_crit
        return 0.0
    else
        return scale * sqrt(kappa - kappa_crit)
    end
end

"""
    simulate_ou_with_resets(N, params, T, dt)

Simulates a toy finite-N system:
    dx_i = (-λ x_i + κ * mean(x)) dt + σ dW_i
and an instantaneous reset to 0 when |x_i| > Θ.

Returns (times, beliefs, means, variances).
"""
function simulate_ou_with_resets(N::Int, params::Dict{Symbol,Float64}, T::Float64, dt::Float64)
    lambda = params[:λ]
    kappa  = params[:κ]
    sigma  = params[:σ]
    theta  = params[:Θ]

    num_steps = floor(Int, T / dt) + 1
    times = collect(0:dt:T)

    beliefs = zeros(Float64, num_steps, N)
    beliefs[1, :] .= randn(N)  # initial condition

    means     = similar(times)
    variances = similar(times)

    means[1]     = mean(beliefs[1, :])
    variances[1] = var(beliefs[1, :])

    sqrt_dt = sqrt(dt)
    for t in 1:(num_steps - 1)
        m = means[t]
        for i in 1:N
            x = beliefs[t, i]
            drift = -lambda * x + kappa * m
            diffusion = sigma * randn() * sqrt_dt

            x_new = x + drift * dt + diffusion

            # instantaneous reset if threshold exceeded
            if abs(x_new) > theta
                x_new = 0.0
            end
            beliefs[t + 1, i] = x_new
        end
        means[t + 1]     = mean(beliefs[t + 1, :])
        variances[t + 1] = var(beliefs[t + 1, :])
    end

    return times, beliefs, means, variances
end

# -----------------------------------------------------------
# Plots
# -----------------------------------------------------------

function plot_bifurcation_diagram(lambda::Float64, sigma::Float64, c0::Float64, nubar::Float64; results_dir::String=".")
    kappa_crit = (sigma^2) / (2 * V_star(lambda, sigma, c0, nubar))
    kappa_values = range(0.0, 1.5 * kappa_crit, length=200)

    stable_variances = [calculate_steady_state_variance(κ, lambda, sigma, c0, nubar) for κ in kappa_values]
    unstable_variances = zeros(length(kappa_values))  # V=0 branch above κ*

    p = plot(kappa_values, stable_variances,
        title = "Supercritical Pitchfork (toy)",
        xlabel = "Peer Influence (κ)",
        ylabel = "Steady-State Variance V*",
        label = "Stable branch",
        lw = 2,
        framestyle = :box,
        size=(900, 550),
    )
    plot!(p, kappa_values, unstable_variances,
        label = "Unstable branch",
        linestyle = :dash,
        lw = 1.5,
    )
    vline!(p, [kappa_crit], label="κ*", linestyle=:dot, color=:black)

    savefig(p, joinpath(results_dir, "bifurcation_diagram.png"))
    return p
end

function plot_phase_portrait(params::Dict{Symbol,Float64}, kappa_to_plot::Float64; results_dir::String=".")
    lambda = params[:λ]
    sigma  = params[:σ]

    p = plot(title="Phase Portrait (κ = $(round(kappa_to_plot, digits=3)))",
             xlabel="Mean μ", ylabel="Variance V",
             framestyle=:box, size=(900, 550))

    # consensus-branch κ* (for coloring switch)
    kappa_crit = (sigma^2) / (2 * V_star(lambda, sigma, 0.0, 0.0))

    # grid
    x_range = range(-2, 2, length=25)
    y_range = range(0,  5, length=25)
    X = repeat(x_range, inner=length(y_range))
    Y = repeat(y_range, outer=length(x_range))

    U = similar(X)              # dμ/dt
    W = similar(Y)              # dV/dt
    dvdV_grid = similar(Y)      # ∂(dV/dt)/∂V
    colors = Vector{Symbol}(undef, length(X))

    for idx in eachindex(X)
        x = X[idx]; y = Y[idx]

        # toy field:
        # dμ/dt = -λ μ (pure mean reversion for μ)
        # dV/dt switches behavior below/above κ*
        U[idx] = -lambda * x
        if kappa_to_plot < kappa_crit
            W[idx] = -y
            dvdV_grid[idx] = -1.0
        else
            W[idx] = -y + y^2
            dvdV_grid[idx] = -1.0 + 2y
        end

        # Diagonal Jacobian eigenvalues: (-λ, dvdV)
        colors[idx] = (max(-lambda, dvdV_grid[idx]) < 0) ? :blue : :red
    end

    stable   = findall(c -> c == :blue, colors)
    unstable = findall(c -> c == :red,  colors)

    quiver!(p, X[stable],   Y[stable],   quiver=(U[stable], W[stable]), color=:blue, label="Stable (Re<0)")
    quiver!(p, X[unstable], Y[unstable], quiver=(U[unstable], W[unstable]), color=:red,  label="Unstable (Re>0)")

    savefig(p, joinpath(results_dir, "ou_phase_portrait.png"))
    return p
end

function plot_time_series(lambda::Float64, sigma::Float64, theta::Float64, N::Int, T::Float64, dt::Float64; results_dir::String=".")
    # κ* from consensus variance (c0=nubar=0 for this diagnostic)
    Vstar0 = V_star(lambda, sigma, 0.0, 0.0)
    kappa_crit = (sigma^2) / (2 * Vstar0)

    # Case A: κ < κ*
    kappa_consensus = 0.5 * kappa_crit
    params_c = Dict(:λ => lambda, :κ => kappa_consensus, :σ => sigma, :Θ => theta)
    times_c, beliefs_c, means_c, vars_c = simulate_ou_with_resets(N, params_c, T, dt)

    # Case B: κ > κ*
    kappa_polar = 1.2 * kappa_crit
    params_p = Dict(:λ => lambda, :κ => kappa_polar, :σ => sigma, :Θ => theta)
    times_p, beliefs_p, means_p, vars_p = simulate_ou_with_resets(N, params_p, T, dt)

    # plot a subset of individuals to keep rendering light
    k = min(N, 50)

    p1 = plot(title="Consensus Regime (κ < κ*)", xlabel="t", ylabel="x_i / μ / V",
              framestyle=:box)
    plot!(p1, times_c, beliefs_c[:, 1:k], label="", color=:lightgray, lw=0.5)
    plot!(p1, times_c, means_c, label="Mean μ(t)", lw=2)
    plot!(p1, times_c, vars_c, label="Variance V(t)", lw=2, yaxis=:right)

    p2 = plot(title="Polarization Regime (κ > κ*)", xlabel="t", ylabel="x_i / μ / V",
              framestyle=:box)
    plot!(p2, times_p, beliefs_p[:, 1:k], label="", color=:lightgray, lw=0.5)
    plot!(p2, times_p, means_p, label="Mean μ(t)", lw=2)
    plot!(p2, times_p, vars_p, label="Variance V(t)", lw=2, yaxis=:right)

    combined = plot(p1, p2, layout=(1, 2), size=(1200, 550))
    savefig(combined, joinpath(results_dir, "time_series_plots.png"))
    return combined
end

function plot_comparative_statics(lambda_fixed::Float64, c0_fixed::Float64, nubar_fixed::Float64; results_dir::String=".")
    # κ* vs σ (recompute V* for each σ)
    sigma_values = range(0.01, 1.0, length=100)
    kappa_star_vs_sigma = [(σ^2) / (2 * V_star(lambda_fixed, σ, c0_fixed, nubar_fixed)) for σ in sigma_values]
    p1 = plot(sigma_values, kappa_star_vs_sigma,
              title="κ* vs Noise σ (with V*(σ))",
              xlabel="σ", ylabel="κ*", lw=2, framestyle=:box, legend=false)

    # κ* vs exogenous V* (for intuition, fixing σ)
    Vstar_values = range(0.01, 0.5, length=100)
    sigma_fixed = 0.5
    kappa_star_vs_V = [(sigma_fixed^2) / (2 * V) for V in Vstar_values]
    p2 = plot(Vstar_values, kappa_star_vs_V,
              title="κ* vs Baseline Dispersion V* (σ fixed)",
              xlabel="V*", ylabel="κ*", lw=2, framestyle=:box, legend=false)

    combined = plot(p1, p2, layout=(1, 2), size=(1200, 550))
    savefig(combined, joinpath(results_dir, "comparative_statics.png"))
    return combined
end

"""
    plot_hopf_bifurcation(lambda, sigma, c0, nubar; results_dir=".")

Toy Hopf-like amplitude diagram with sqrt scaling above κ*.
"""
function plot_hopf_bifurcation(lambda::Float64, sigma::Float64, c0::Float64, nubar::Float64; results_dir::String=".")
    Vstar = V_star(lambda, sigma, c0, nubar)
    kappa_crit = (sigma^2) / (2 * Vstar)
    kappa_values = range(0.0, 1.5 * kappa_crit, length=200)
    amplitudes = [κ <= kappa_crit ? 0.0 : sqrt(κ - kappa_crit) for κ in kappa_values]

    p = plot(kappa_values, amplitudes,
             title="Hopf-like Amplitude (toy)",
             xlabel="Peer Influence (κ)",
             ylabel="Limit Cycle Amplitude",
             lw=2, framestyle=:box, label="Amplitude", size=(900, 550))
    vline!(p, [kappa_crit], label="κ*", linestyle=:dot, color=:black)

    savefig(p, joinpath(results_dir, "hopf_bifurcation.png"))
    return p
end

# -----------------------------------------------------------
# Orchestration
# -----------------------------------------------------------

function generate_ou_with_resets_plots(results_dir::String)
    mkpath(results_dir)
    Random.seed!(42)  # reproducibility

    # Baseline parameters (toy)
    lambda = 0.5
    sigma  = 0.1
    theta  = 1.0
    N      = 200
    T      = 100.0
    dt     = 0.1
    c0     = 0.05    # contraction factor (only used in V_star)
    nubar  = 0.15    # mean jump rate proxy (only used in V_star)

    # Diagnostics for κ*
    Vstar_val       = V_star(lambda, sigma, c0, nubar)
    kappa_crit_val  = (sigma^2) / (2 * Vstar_val)

    # 1) Bifurcation-style diagram
    plot_bifurcation_diagram(lambda, sigma, c0, nubar; results_dir=results_dir)

    # 2) Phase portrait at κ = 0.8 κ*
    plot_phase_portrait(Dict(:σ => sigma, :λ => lambda), 0.8 * kappa_crit_val; results_dir=results_dir)

    # 3) Time series in consensus vs polarization regimes
    plot_time_series(lambda, sigma, theta, N, T, dt; results_dir=results_dir)

    # 4) Comparative statics for κ*
    plot_comparative_statics(lambda, c0, nubar; results_dir=results_dir)

    # 5) Hopf-like amplitude diagram
    plot_hopf_bifurcation(lambda, sigma, c0, nubar; results_dir=results_dir)

    println("✅ Plots saved in: $(abspath(results_dir))")
end

function main()
    results_dir = get(ENV, "RESULTS_DIR", "results")
    generate_ou_with_resets_plots(results_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
