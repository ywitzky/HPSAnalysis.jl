using HDF5, StatsBase

Sim = HPSAnalysis.SimData()
Sim.SlabAxis=3
Sim.BoxLength = [10,10,10]
Sim.BasePath = TestPath * "SurfaceTensionTest/"

testname ="testfile.h5"
mkpath(TestPath * "SurfaceTensionTest/")
filename=TestPath * "SurfaceTensionTest/$(testname)"
testfile=fid = h5open(filename, "w")

key = "hoomd-data/md/compute/ThermodynamicQuantities/pressure_tensor"
N = 10_202
pressure_tensor = rand(6,N)
testfile[key] = pressure_tensor

timestep = (1:N) .* 100
key = "hoomd-data/Simulation/timestep"
testfile[key] = collect(timestep)


start = 20
NFrames = 1000
Sim.TrajWriteOutFreq =NFrames
τ=2
NSub=5
γ, Δγ = HPSAnalysis.computeSurfaceTension(Sim,Int32(start);filename=testname, tau=Int32(τ), NSub=Int32(NSub))

gamma_vec = Sim.BoxLength[Sim.SlabAxis]/2.0 *(pressure_tensor[6,200:τ:end].-(pressure_tensor[1,200:τ:end].+pressure_tensor[4,200:τ:end])./2.0)
gamma_mean = mean(gamma_vec)

gamma_sub = [mean(gamma_vec[a:e]) for (a,e) in [(1, 1000),(1001, 2000),(2001, 3000),(3001, 4000),(4001, 5000)]]
gamma_error = sqrt(sum((gamma_sub .- gamma_mean).^2)/5)

@testset "Surface Tension" begin
    @test γ  ≈ gamma_mean
    @test Δγ ≈ gamma_error
end