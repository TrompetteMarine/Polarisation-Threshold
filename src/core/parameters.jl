using Parameters

@with_kw struct ModelParams
    λ::Float64 = 1.0
    σ::Float64 = 0.5
    κ::Float64 = 0.1
    Θ::Float64 = 1.0
    ϖ::Float64 = 0.1
    N::Int = 1000
    T::Float64 = 100.0
    dt::Float64 = 0.01
    seed::Int = 1234
end
