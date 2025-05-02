using GSDFormat

SetupTestPath="$(TestPath)/Setup/"
rm(SetupTestPath; force=true, recursive=true)
mkpath(SetupTestPath)
pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
sim = pyimport("Submit_HOOMD")
###############


pH=8.0
NChains = 2
Names=["A","B"]
Sequences=["EEEEE", "DDDDD"] 
BoxLengthShort=200.0
BoxLengthLong=2000.0
BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]
Info="Test Sims"

pos = zeros(Float32, 10, 3)

pos[1:5, 2] .= collect(1:5)*3.8
pos[6:10, 2] .= collect(1:5)*3.8.+25

println(SetupTestPath)
cd(SetupTestPath)
HPSAnalysis.Setup.writeStartConfiguration("./Test_slab","./Test_Start_slab.txt", Info, Sequences, BoxSize , 10_000, HOOMD=true, ; SimulationType="Calvados2" , Temperature=300,  InitStyle="Pos", Pos=pos , pH=pH, SaltConcentration=0.2, Device="CPU", WriteOutFreq=1_000)

sim.run(SetupTestPath)
sim.restart(SetupTestPath; ExtendedSteps=20_000)
sim.restart(SetupTestPath; ExtendedSteps=30_000)
sim.restart(SetupTestPath; ExtendedSteps=40_000)
sim.restart(SetupTestPath; ExtendedSteps=50_000)

traj_1 = GSDFormat.open("$(SetupTestPath)/traj.gsd")
traj_2 = GSDFormat.open("$(SetupTestPath)/traj_1.gsd")
traj_3 = GSDFormat.open("$(SetupTestPath)/traj_2.gsd")
traj_4 = GSDFormat.open("$(SetupTestPath)/traj_3.gsd")
traj_5 = GSDFormat.open("$(SetupTestPath)/traj_4.gsd")

data = HPSAnalysis.initData(SetupTestPath; Reparse=true, LoadAll=true, HOOMD=true, ReadBig=true)

@testset "Restarts" begin
    @test data.NSteps == 50 

    for (traj, range) in zip([traj_1, traj_2, traj_3, traj_4, traj_5], [1:10, 11:20, 21:30,31:40,41:50]) 
        for (step, frame) in zip(range,traj)
            @test all(frame.particles.position[:, 1].*10.0 .≈ data.x[:, step])
            @test all(frame.particles.position[:, 2].*10.0 .≈ data.y[:, step])
            @test all(frame.particles.position[:, 3].*10.0 .≈ data.z[:, step])
        end
    end
end