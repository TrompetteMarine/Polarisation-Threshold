using DifferentialEquations

function ou_drift!(du, u, p, t)
    λ, κ, G = p.λ, p.κ, p.G
    du[1] = -λ * u[1] + κ * G
    du[2] = 0.0
end

function ou_diffusion!(du, u, p, t)
    σ = p.σ
    du[1, 1] = σ
    du[2, 1] = 0.0
end
