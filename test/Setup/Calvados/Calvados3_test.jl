SimName="test"

import HPSAnalysis.BioData as BioData
using GSDFormat 

if isdir("$SetupTestPath/HOOMD_Setup/")
    rm("$SetupTestPath/HOOMD_Setup/"; force=true, recursive=true)
end
mkpath("$SetupTestPath/HOOMD_Setup/")


Seq=["DEGHKDEGHK"]
HPSAnalysis.Setup.WriteHOOMDSequences("$SetupTestPath/HOOMD_Setup/Sequences.txt", Seq)

N=10
NChains=1
InputBonds=[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9]]
InputAngles=[]
InputDihedrals=[]
AtomTypes= Set(join(Seq))
AaToId = Dict{Char,Int32}()
for (index, value) in enumerate(AtomTypes)
    AaToId[value]=index
end


position=[0.0 0.0 1.0 0.5 0.0 0.0 -0.5 -1.0 0.0 0.0;;;-4.0 -3.0 -2.0 -1.0 0.0 0.0 1.0 2.0 3.0 4.0;;;0.0 0.0 -1.0 -0.5 -1.0 1.0 0.5 1.0 0.0 0.0]


coor=fill(0.0,1,maximum(length(seq) for seq in Seq),3)



HPSAnalysis.Setup.WriteHOOMDParticlesInput("$SetupTestPath/HOOMD_Setup/Particles.txt",position,BioData.OneToHPSCharge,AaToId,Seq,BioData.AaToWeight,BioData.OneToHPSCalvadosSigma,coor)
HPSAnalysis.Setup.WriteDictionaries("$SetupTestPath/HOOMD_Setup/Dictionaries.txt", BioData.OneToHPSCharge, AaToId,BioData.AaToWeight, BioData.OneToHPSCalvadosSigma, BioData.OneToCalvados2Lambda)
HPSAnalysis.Setup.WriteParams("$SetupTestPath/HOOMD_Setup/Params.txt",SimName,300, 10, 1, 0.01, Array([10,101,10]), rand(1:65535), UseAngles=false,domain=Array([[3,8]]),Device="CPU", UseCharge=false, Create_Start_Config=true,SimType="Calvados3")
HPSAnalysis.Setup.WriteDihedrals("$SetupTestPath/HOOMD_Setup/DihedralMap.txt",[],0)

#@error "ENM are not properly tested!"
harmonic = Dict{String, Dict{Symbol, Float64}}()
r0 = 5
harmonic["O-O"] = Dict(:r => r0, :k => 700)
harmonic["BB_1"] = Dict(:r => r0, :k => 700)
harmonic["ENM_2"] = Dict(:r => r0, :k => 1800)
harmonic["ENM_3"] = Dict(:r => r0, :k => 1800)
harmonic["ENM_4"] = Dict(:r => r0, :k => 1800)
HPSAnalysis.Setup.WriteENM_HOOMD_Indices("$SetupTestPath/HOOMD_Setup/ENM_indices.txt", (12, ["O-O", "BB_1", "ENM_2", "ENM_3", "ENM_4"], [1, 1, 1, 2, 1, 1, 1, 1, 1, 3, 4, 5], [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (3, 5), (3, 6), (4, 6)], harmonic))

sim.run(SetupTestPath)
data=GSDFormat.open("$(SetupTestPath)$(SimName)_300.0_Start_slab.gsd","r")


frame = data[1]
particle_N_test = frame.particles.N
particle_position_test = frame.particles.position
particle_types_test = frame.particles.types
particle_typeid_test = frame.particles.typeid
particle_image_test = frame.particles.image

bond_N_test=frame.bonds.N
bond_types_test=frame.bonds.types
bond_typid_test=frame.bonds.typeid
bond_group_test=frame.bonds.group

typesid=UInt32[1, 3, 4, 2, 0, 1, 3, 4, 2, 0]
coor=Float32[0.0 -0.4 0.0; 0.0 -0.3 0.0; 0.1 -0.2 -0.1; 0.05 -0.1 -0.05; 0.0 0.0 -0.1; 0.0 0.0 0.1; -0.05 0.1 0.05; -0.1 0.2 0.1; 0.0 0.3 0.0; 0.0 0.4 0.0]
bond_group=Int32[0 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 3 5; 3 6; 4 6]
image=Int32[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
types=["K", "D", "H", "E", "G"]

bondid=UInt32[0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 3, 4]
bondtypes=["O-O\0\0", "BB_1\0", "ENM_2", "ENM_3", "ENM_4"]

@testset "Calvados3" begin
    @test particle_N_test==10
    @test particle_position_test==coor
    @test particle_types_test==types
    @test particle_typeid_test==typesid
    @test particle_image_test==image
    
    @test bond_N_test==UInt32(12)
    @test bond_types_test==bondtypes
    @test bond_typid_test==bondid
    @test bond_group_test==bond_group
end

close(data)
