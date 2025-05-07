pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
sim = pyimport("Submit_HOOMD")
mkpath("$SetupTestPath/HOOMD_Setup/")


@testset "Calvados Tests" begin
    #include("./Calvados2_test.jl")
    #include("./Calvados3_test.jl")
    include("RS_Prot.jl")
end
