filename_test="$(SetupTestPath)/GSD_write_test.gsd"
filename="$(SetupTestPath)/GSD_write.gsd"

@testset "writeGSDStartFile" begin
    sequences=["MNAM"]
    N=4
    set=Set(join(sequences))
    AAToID=Dict{Char,Int}()
    max_seq_length = 4
    BoxSize=Vector{Float64}([20.0,30.0,20.0,0.0,0.0,0.0])
    coor = HPSAnalysis.Setup.createStartingPosition(sequences,BoxSize)
    AltBox = [BoxSize[2]-BoxSize[1], BoxSize[4]-BoxSize[3], BoxSize[6]-BoxSize[5]]
    coor = HPSAnalysis.Setup.correctPositionInBounds(coor, AltBox, sequences)
    coor = HPSAnalysis.Setup.correctPositionInBounds(coor, AltBox, sequences)
    coor = HPSAnalysis.Setup.correctPositionInBounds(coor, AltBox, sequences)
    for (num,atom) in enumerate(set)
        AAToID[atom]=num
    end
    image = HPSAnalysis.Setup.getImageCopyNumber(coor, AltBox, sequences)

    snapshot=GSDFormat.Frame()
    snapshot.configuration.step=1
    snapshot.configuration.dimensions=3
    snapshot.configuration.box=[20,30,20,0,0,0]./10.0
    snapshot.particles.N=N
    snapshot.particles.position=reshape(permutedims(coor, (2,1,3)), (size(coor, 1)*size(coor, 2), 3))./10.0
    IDToAA=Dict(value=>key for (key,value) in AAToID)
    snapshot.particles.types =  [string(IDToAA[Id]) for Id in  sort(collect(values(AAToID))) ]
    snapshot.particles.typeid = [Int32(AAToID[AA])-1 for AA in join(sequences)]
    snapshot.particles.image = reshape(permutedims(image, (2,1,3)), (size(image, 1)*size(image, 2), 3))
    mass_charge=Vector{Float64}([1.0,2.0])
    snapshot.particles.mass = mass_charge
    snapshot.particles.charge = mass_charge
    snapshot.particles.diameter = [Float32(AAToID[AA])/10.0  for AA in join(sequences)]

    snapshot.bonds.N=N-1
    snapshot.bonds.types = ["O-O"]
    snapshot.bonds.typeid = zeros(Int32, N-1)
    snapshot.bonds.group = HPSAnalysis.Setup.getBonds(sequences, M=2)

    UseAngles=false
    DiMap=Dict()
    DiList=Matrix{Int}(undef,4,4)
    Ditypes=[]
    if UseAngles
        # Create Angles
        snapshot.angles.N = N-1
        snapshot.angles.types = ["O-O-O"]
        snapshot.angles.typeid = zeros(Int32, N-2)
        snapshot.angles.group = HPSAnalysis.Setup.getBonds(sequences, M=3)

        # Create Dihedrals
        snapshot.dihedrals.N = N-3

        snapshot.dihedrals.types = ["$(ids[1])-$(ids[2])-$(ids[3])-$(ids[4])" for ids in collect(keys(DihedralMap))]#string.(collect(values(DihedralMap)))
        snapshot.dihedrals.typeid = [DihedralMap[DihedralList[key,:]]-1 for key in axes(DihedralList,1)] ### convert to python numbering
        snapshot.dihedrals.group = HPSAnalysis.Setup.getBonds(Sequences, M=4)
    end

    file = GSDFormat.open(filename_test, 'w')
    GSDFormat.append(file, snapshot)
    GSDFormat.close(file)

    HPSAnalysis.Setup.writeGSDStartFile(filename,N,N-1,N-2,N-3,BoxSize,coor,AAToID,sequences,image,mass_charge,mass_charge,DiMap,DiList,AAToID,false)

    @test files_are_equal(filename_test,filename)
end
rm(filename_test)
rm(filename)