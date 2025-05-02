SetupTestPath="$(TestPath)/Setup/"
mkpath(SetupTestPath)

function files_are_equal(file1::String, file2::String)::Bool
    return read(file1, String) == read(file2, String)
end

@testset "Setup Tests" begin
    include("./Test_setup.jl")
    include("./HOOMD_Setup_test.jl")
    include("./GSD_test.jl")
    include("./Restart_test.jl")
end
