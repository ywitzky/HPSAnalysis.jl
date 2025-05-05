using HPSAnalysis
PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])

EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)
#=
file="$(PkgSourcePath)../data/EnvironmentPath.txt" 
if isfile(file)
    f = open(file,"r")
    EnvironmentPath = readline(f)
end
if EnvironmentPath==""
    @warn("Python Environment has not been specified. All tests regarding Polyply and HOOMD scripts will be omitted.")
end=#

ENV["PYCALL_JL_RUNTIME_PYTHON"]="$(EnvironmentPath)/bin/python"

using PyCall, Test, Scratch, Aqua
pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")

TestPath = Scratch.get_scratch!(HPSAnalysis, "test") 

#@testset "Aqua" begin
#    Aqua.test_all(HPSAnalysis; deps_compat=(ignore=[:Printf, :Mmap, :Libdl, :LinearAlgebra, :Statistics,:Test], ), project_extras=false, )
    ### ignore standard libraries, not sure how to deal with them im PackageCompatUI/add compats manually
    ### Test fails in project extras since they dont get excluded their normaly
#end


#include("./Analysis/Analysis_test.jl")
include("./Setup/Setup_test.jl")
include("./Calvados/C_test.jl")
