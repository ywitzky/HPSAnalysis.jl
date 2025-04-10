SimName="test"

Seq=["DEGHKDEGHK"]
HPSAnalysis.Setup.WriteHOOMDSequences("$TestPath/HOOMD_Setup/Sequences.txt", Seq)

N=10
NChains=1
InputBonds=[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]]
InputAngles=[]
InputDihedrals=[]
AtomTypes= Set(join(Seq))
AaToId = Dict{Char,Int32}()
for (index, value) in enumerate(AtomTypes)
    AaToId[value]=index
end


position=[0.0 0.0 1.0 0.5 0.0 0.0 -0.5 -1.0 0.0 0.0;;;-4.0 -3.0 -2.0 -1.0 0.0 0.0 1.0 2.0 3.0 4.0;;;0.0 0.0 -1.0 -0.5 -1.0 1.0 0.5 1.0 0.0 0.0]
#Array([[[0.0,-4.0,0.0],[0.0,-3.0,0.0],[1.0,-2.0,-1.0],[0.5,-1.0,-0.5],[0.0,0.0,-1.0],[0.0,0.0,1.0],[-0.5,1.0,0.5],[-1.0,-2.0,1.0],[0.0,3.0,0.0],[0.0,4.0,0.0]]])

coor=fill(0.0,1,maximum(length(seq) for seq in Seq),3)


HPSAnalysis.Setup.WriteHOOMDParticlesInput("$TestPath/HOOMD_Setup/Particles.txt",position,BioData.OneToHPSCharge,AaToId,Seq,BioData.AaToWeight,BioData.OneToHPSCalvadosSigma,coor)
HPSAnalysis.Setup.WriteDictionaries("$TestPath/HOOMD_Setup/Dictionaries.txt", BioData.OneToHPSCharge, AaToId,BioData.AaToWeight, BioData.OneToHPSCalvadosSigma, BioData.OneToCalvados2Lambda)
HPSAnalysis.Setup.WriteParams("$TestPath/HOOMD_Setup/Params.txt",SimName,300, 10, 1, 0.01, Array([10,101,10]), rand(1:65535), UseAngles=false,domain=Array([[3,8]]),Device="CPU", UseCharge=false, Create_Start_Config=true,SimType="Calvados3")
HPSAnalysis.Setup.WriteDihedrals("$TestPath/HOOMD_Setup/DihedralMap.txt",[],0)

sim.run("$TestPath")
data=GSDFormat.open("$(TestPath)$(SimName)_StartConfiguration.gsd","r")

for (i,frame) in enumerate(data)
    global particle_N_test = frame.particles.N
    global particle_position_test = frame.particles.position
    global particle_types_test = frame.particles.types
    global particle_typeid_test = frame.particles.typeid
    global particle_image_test = frame.particles.image

    global bond_N_test=frame.bonds.N
    global bond_types_test=frame.bonds.types
    global bond_typid_test=frame.bonds.typeid
    global bond_group_test=frame.bonds.group
end

typesid=UInt32[1, 3, 4, 2, 0, 1, 3, 4, 2, 0]
coor=Float32[0.0 -0.4 0.0; 0.0 -0.3 0.0; 0.1 -0.2 -0.1; 0.05 -0.1 -0.05; 0.0 0.0 -0.1; 0.0 0.0 0.1; -0.05 0.1 0.05; -0.1 0.2 0.1; 0.0 0.3 0.0; 0.0 0.4 0.0]
bond_group=Int32[0 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 3 5; 3 6; 3 7; 4 6; 4 7; 5 7]
image=Int32[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
types=["K", "D", "H", "E", "G"]

bondid=UInt32[0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
bondtypes=["O-O\0\0\0\0", "O-O_D1\0", "O-O_D2\0", "O-O_D3\0", "O-O_D4\0", "O-O_D5\0", "O-O_D6\0", "O-O_D7\0", "bond_8\0", "bond_9\0", "bond_10", "bond_11", "bond_12", "bond_13"]
#["O-O", "O-O_D1", "O-O_D2", "O-O_D3", "O-O_D4", "O-O_D5", "O-O_D6", "O-O_D7", "bond_8", "bond_9", "bond_10", "bond_11", "bond_12", "bond_13"]

@testset "Calvados3" begin
    @test particle_N_test==10
    @test particle_position_test==coor
    @test particle_types_test==types
    @test particle_typeid_test==typesid
    @test particle_image_test==image
    
    @test bond_N_test==UInt32(15)
    @test bond_types_test==bondtypes
    @test bond_typid_test==bondid
    @test bond_group_test==bond_group
end
