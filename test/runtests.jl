using HPSAnalysis, Aqua, Test, Scratch

PkgPath=pathof(HPSAnalysis)
TestPath = Scratch.get_scratch!(HPSAnalysis, "test") #"$PkgPath/test/TemporaryFiles/"

@testset "Aqua" begin
    Aqua.test_all(HPSAnalysis; deps_compat=(ignore=[:Printf, :Mmap, :Libdl, :LinearAlgebra, :Statistics,:Test], ), project_extras=false, )
    ### ignore standard libraries, not sure how to deal with them im PackageCompatUI/add compats manually
    ### Test fails in project extras since they dont get excluded their normaly
end

include("./Analysis/Analysis_test.jl")
include("./Setup/Setup_test.jl")