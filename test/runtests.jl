using HPSAnalysis
PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
EnvironmentPath="/localscratch/Programs/hoomd/.cpu2/bin/python3"
ENV["PYCALL_JL_RUNTIME_PYTHON"]=EnvironmentPath

using PyCall, Test, Scratch#, Aqua
pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")

TestPath = Scratch.get_scratch!(HPSAnalysis, "test") #"$PkgPath/test/TemporaryFiles/"
#@testset "Aqua" begin
#    Aqua.test_all(HPSAnalysis; deps_compat=(ignore=[:Printf, :Mmap, :Libdl, :LinearAlgebra, :Statistics,:Test], ), project_extras=false, )
    ### ignore standard libraries, not sure how to deal with them im PackageCompatUI/add compats manually
    ### Test fails in project extras since they dont get excluded their normaly
#end


#include("./Analysis/Analysis_test.jl")
include("./Setup/Setup_test.jl")
include("./Calvados/C_test.jl")
