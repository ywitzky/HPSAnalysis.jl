SetupTestPath="$(TestPath)/Setup/"
mkpath(SetupTestPath)

pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
sim = pyimport("Submit_HOOMD")
mkpath("$SetupTestPath/HOOMD_Setup/")

function files_are_equal(file1::String, file2::String)::Bool
    return read(file1, String) == read(file2, String)
end#=
include("./Calvados/Calvados3_ENM_test.jl")

@testset "Setup Tests" begin
    include("./HOOMD_Setup_test.jl")
    include("./GSD_test.jl")
end
=#
@testset "Calvados Tests" begin
    #include("./Calvados/Calvados_Helper_test.jl")
    #include("./Calvados/Calvados2_test.jl")
    #include("./Calvados/Calvados3_test.jl")
    include("./Calvados/Implementation/Calvados3_implementation_test.jl")
    include("./Calvados/Implementation/Calvados3_implementation_analysis.jl")

    #include("RS_Prot.jl") not yet fully implemented
end
