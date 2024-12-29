module LammpsAnalysis

### helping stuff which probably should be seperate
include("./BioData.jl")
include("./ProteinSequences.jl")
#import .BioData
#import .ProteinSequences
include("./Xtc.jl") 
include("./Setup.jl")

 
#include("./StructDispatch.jl") ###TODO: Do proper typing and get rid of this at some point
include("./Data.jl")
include("./Polyply.jl")
include("./IO.jl")
include("./Analysis.jl")
include("./Unify.jl")
include("./Plot.jl")
include("./HREMD.jl")
end # module LammpsAnalysis
