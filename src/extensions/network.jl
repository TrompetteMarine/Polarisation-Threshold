function create_adjacency_matrix(N::Int, p::Float64)
    A = rand(N, N) .< p
    A = A .* (1 .- I(N))  # No self-loops
    return A
end

function compute_mean_field(agents::Vector{Agent}, A::Matrix{Float64}, i::Int)
    return sum(A[i, j] * (agents[j].x - agents[j].r) for j in 1:length(agents)) / sum(A[i, j] for j in 1:length(agents))
end
