module HPSAnalysis

include("./IO/PythonEnv.jl")
PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
PkgPath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-2])
EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)


### helping stuff which probably should be seperate
include("./BioData.jl")
include("./ProteinSequences.jl")
include("./Data.jl")

#import .BioData
#import .ProteinSequences
include("./Xtc.jl") 
include("./Setup.jl")

 
#include("./StructDispatch.jl") ###TODO: Do proper typing and get rid of this at some point
include("./Polyply.jl")
include("./IO.jl")
include("./Analysis.jl")
include("./Plot.jl")
end # module HPSAnalysis
