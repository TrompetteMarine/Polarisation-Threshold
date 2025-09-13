module ModelCore

include("parameters.jl")
include("agent.jl")
include("dynamics.jl")

export ModelParams, Agent, initialize_agents, ou_drift!, ou_diffusion!

end # module ModelCore
