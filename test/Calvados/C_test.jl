pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
sim = pyimport("Submit_HOOMD")

@testset "Calvados Tests" begin
    include("./Calvados2_test.jl")
    include("./Calvados3_test.jl")
end
