using GSDFormat, HDF5

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

ENM = (0,[],[],[], Dict())
cd(SetupTestPath)
HPSAnalysis.Setup.writeStartConfiguration("./Test_slab","./Test_Start_slab.txt", Info, Sequences, BoxSize , 10_000, HOOMD=true, ; SimulationType="Calvados2" , Temperature=300,  InitStyle="Pos", Pos=pos , pH=pH, SaltConcentration=0.2, Device="CPU", WriteOutFreq=1_000, ENM=ENM)

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

@testset "Restarts" begin
    ### test simulation extension by restarts
    data = HPSAnalysis.initData(SetupTestPath; Reparse=true, LoadAll=true, HOOMD=true, ReadBig=true)
    @test data.NSteps == 50 

    #=
    ### test if data is correctly read by initData
    alltrue = true ### reduce the number of tests since basically one test
    for (traj, range) in zip([traj_1, traj_2, traj_3, traj_4, traj_5], [1:10, 11:20, 21:30,31:40,41:50]) 
        for (step, frame) in zip(range,traj)
            alltrue &= all(frame.particles.position[:, 1].*10.0 .≈ data.x[:, step])
            alltrue &= all(frame.particles.position[:, 2].*10.0 .≈ data.y[:, step])
            alltrue &= all(frame.particles.position[:, 3].*10.0 .≈ data.z[:, step])
        end
    end
    @test alltrue
    =#

    println(SetupTestPath)
    file = h5open("$(SetupTestPath)/pressure.h5")
    steps = file["hoomd-data"]["Simulation"]["timestep"][:]
    @test length(steps) == 500
    expected_steps = vcat(collect(8100:100:18000), collect(100:100:10_000),collect(100:100:10_000) ,collect(100:100:10_000) ,collect(100:100:10_000) )
    @test all(steps .== expected_steps)
    close(file)

    ### restart "crashed" simulation
    HPSAnalysis.Setup.WriteParams("$(SetupTestPath)/HOOMD_Setup/Params.txt", "./Test_Start_slab.txt", 300, 60_000, 1_000, 0.01, [20.0,200.0,20.0], 11111; UseAngles=false, Device="CPU")

    sim.restart(SetupTestPath; ExtendedSteps=0) ### runs a new traj for additional 10_000 steps up to 60_000
    file = h5open("$(SetupTestPath)/pressure.h5")
    steps = file["hoomd-data"]["Simulation"]["timestep"][:]
    traj_6 = GSDFormat.open("$(SetupTestPath)/traj_5.gsd")

    @test length(steps) == 600
    @test all(steps .== vcat(expected_steps, collect(100:100:10000)))
    data = HPSAnalysis.initData(SetupTestPath; Reparse=true, LoadAll=true, HOOMD=true, ReadBig=true)

    ### test if data is correctly read by initData
    alltrue = true ### reduce the number of tests since basically one test
    for (traj, range) in zip([traj_1, traj_2, traj_3, traj_4, traj_5,traj_6], [1:10, 11:20, 21:30,31:40,41:50, 51:60]) 
        for (step, frame) in zip(range,traj)
            alltrue &= all(frame.particles.position[:, 1].*10.0 .≈ data.x[:, step])
            alltrue &= all(frame.particles.position[:, 2].*10.0 .≈ data.y[:, step])
            alltrue &= all(frame.particles.position[:, 3].*10.0 .≈ data.z[:, step])
        end
    end
    @test alltrue
end


