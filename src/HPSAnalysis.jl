module HPSAnalysis

### helping stuff which probably should be seperate
include("../data/BioData.jl")
include("../data/ProteinSequences.jl")

include("./Setup.jl")

 
include("./Data.jl")
include("./Polyply.jl")
include("./IO.jl")
include("./Analysis.jl")
include("./Unify.jl")
include("./Plot.jl")
include("./HREMD.jl")
end # module HPSAnalysis
