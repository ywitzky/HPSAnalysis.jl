SetupTestPath="$(TestPath)/Setup/"
rm(SetupTestPath)
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
