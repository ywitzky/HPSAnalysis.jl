
@testset "Analysis Tests" begin
    include("./Slab_Analysis_test.jl")
    include("./Clustering_test.jl")
    include("./Slab_Cumulant_test.jl")
    include("./IntraChainScaling_test.jl")
    include("./IntraChainContactMatrix_test.jl")
end