module PolarisationThreshold

include("Parameters.jl")
include("ABM/ABM.jl")
include("MeanField/MeanField.jl")
include("Analysis/Analysis.jl")
include("Visualization/Visualization.jl")

using .Parameters
using .ABM
using .MeanField
using .Analysis
using .Visualization

export Parameters, ABM, MeanField, Analysis, Visualization
export ModelParams, num_steps, as_namedtuple, update

const ModelParams = Parameters.ModelParams
const num_steps = Parameters.num_steps
const as_namedtuple = Parameters.as_namedtuple
const update = Parameters.update

end # module PolarisationThreshold
