# Polarisation Threshold Simulation Toolkit

The PolarisationThreshold package implements stochastic and deterministic models of opinion formation using the Ornstein–Uhlenbeck-with-resets mechanism. The code base is organised as a Julia package so that you can run agent-based simulations, analyse mean-field dynamics, generate publication-ready figures, and collect summary metrics from a consistent set of parameters.

## Repository Layout

The `src/PolarisationThreshold.jl` entry point re-exports five cohesive submodules:

- **Parameters** – provides the `ModelParams` struct and helpers such as `num_steps`, `as_namedtuple`, and `update` for manipulating parameter sets.
- **ABM** – defines the `OpinionAgent` type, `initialize_model`, and the `run_abm` orchestration helper built on top of Agents.jl.
- **MeanField** – exposes `mean_field_system`, `reduced_reset_map`, `orbit_data`, and `bifurcation_data` utilities backed by DynamicalSystems.jl.
- **Analysis** – contains metrics that compare stochastic and deterministic observables (e.g., `kappa_star_from_abm`, `compare_thresholds`).
- **Visualization** – wraps Plots.jl recipes for ABM time series and deterministic orbit/bifurcation outputs.

Each submodule can be imported directly via `using PolarisationThreshold`.

## Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/Polarisation-Threshold.git
   cd Polarisation-Threshold
   ```

2. **Instantiate the Julia project**
   ```bash
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```

This installs Agents.jl, AgentsPlots.jl, DynamicalSystems.jl, ChaosTools.jl (used via DynamicalSystems), Plots.jl, JLD2, and the rest of the declared dependencies.

## Working with Parameters

The `ModelParams` struct centralises all tunable values:

```julia
using PolarisationThreshold
params = ModelParams(λ = 1.0, σ = 0.5, κ = 0.2, Θ = 1.0, ϖ = 0.2, N = 500, T = 200.0, dt = 0.01, seed = 1234)
```

Use `update(params; kwargs...)` to create modified copies and `num_steps(params)` to convert `(T, dt)` into a discrete step count. Converting to a `NamedTuple` (`as_namedtuple(params)`) is convenient when logging configurations.

## Running the Agent-Based Model

```julia
using PolarisationThreshold
params = ModelParams(N = 250, κ = 0.18)
result = ABM.run_abm(params; steps = 10_000)
```

`run_abm` initialises an `Agents.AgentBasedModel`, advances it for the requested number of steps (defaults to `num_steps(params)`), and returns a named tuple containing:

- `model` – the evolved `AgentBasedModel` instance.
- `agent_data` – raw agent-level data (when `store = true`).
- `model_data` – the backend-specific collector payload.
- `data` – a normalised `Dict{Symbol, Vector}` of observables (e.g. `:G`, `:mean_x`, `:variance_u`).
- `observables` – the mapping that was passed to the collector.

Custom observables can be supplied via a dictionary or iterable of `(name, function)` pairs. The helpers in `Visualization` accept the `data` dictionary and render plots such as

```julia
Visualization.plot_model_observables(result.data)
Visualization.plot_mean_field_series(result.data)
Visualization.plot_variance_series(result.data)
```

All figures are written to the `results/` directory and the returned `Plots.Plot` objects can be further customised.

## Deterministic Mean-Field Workflows

Construct the continuous mean-field system or reduced reset map to analyse deterministic limits:

```julia
params = ModelParams(κ = 0.25)
system = MeanField.mean_field_system(params)
reduced = MeanField.reduced_reset_map(params; initial_state = [0.1])
```

Generate trajectories with `orbit_data(system; steps = 2_000, discard = 200)` and sweep parameter ranges using

```julia
κ_values = 0.0:0.01:0.5
factory = κ -> MeanField.mean_field_system(update(params; κ = κ); reduced = true)
bif = MeanField.bifurcation_data(factory, κ_values; observable = state -> state[1], steps = 1_200, discard = 200)
```

Visualisations are available via `Visualization.plot_orbit_trajectory(orbit)` and `Visualization.plot_bifurcation_diagram(bif)`, again saving to `results/`.

## Analysis Utilities

`src/Analysis/metrics.jl` provides quick diagnostic routines for comparing stochastic and deterministic behaviour:

- `variance_stream(model_data)` – extracts the variance observable from ABM results.
- `kappa_star_from_abm(model_data, params; window=10)` – estimates the critical coupling from the tail variance.
- `stability_fraction(orbit, params)` – measures the time spent inside the stability threshold in deterministic orbits.
- `kappa_star_from_orbit(orbit, params)` – deterministic analogue of the κ* estimate.
- `compare_thresholds(abm_data, orbit, params)` – returns a named tuple with stochastic vs deterministic κ* and their difference.

These functions operate directly on the outputs of `ABM.run_abm` and `MeanField.orbit_data`.

## Command-Line Workflows

Two scripts in `scripts/` orchestrate common experiments:

- `scripts/run_abm.jl` instantiates the project, executes the ABM with default parameters, writes plots to `results/`, and stores the processed observable dictionary in `results/abm_results.jld2`.
- `scripts/run_dynamics.jl` generates a deterministic orbit, sweeps κ over a configurable range to build a bifurcation diagram, saves plots, and persists the raw orbit plus bifurcation data to `results/mean_field_results.jld2`.

Run them with the Julia executable on your system:

```bash
julia --project scripts/run_abm.jl
julia --project scripts/run_dynamics.jl
```

Both scripts activate the package environment automatically and extend `LOAD_PATH` so that `PolarisationThreshold` resolves correctly when executed outside the package context.

## Output Artefacts

All plotting functions create the `results/` directory if necessary and save PNG files there. The orchestration scripts also persist raw data (JLD2 format) that can be reloaded for documentation builds or post-processing.

## Testing

Execute the package test suite with

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

The tests exercise ABM setup, mean-field constructors, and analysis routines to ensure the documented workflows remain operational.

