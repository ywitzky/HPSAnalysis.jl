filename_test="$(SetupTestPath)/GSD_write_manual.gsd"
filename="$(SetupTestPath)/GSD_write_setup.gsd"

using GSDFormat

macro namedtest(name, test)
    esc(:(@testset $name begin @test $test end))
end

@testset "writeGSDStartFile" begin
    #Not testing UseAngles=true, because it is an artefact of a previous version.
    UseAngles=false

    sequences=["MNAM"]
    N=4
    set=Set(join(sequences))
    AAToID=Dict{Char,Int}()
    max_seq_length = 4
    BoxSize=Vector{Float64}([20.0,30.0,20.0,0.0,0.0,0.0])
    coor_notreshaped=[5.0 5.0 5.0 5.0;;; 62.4 66.2 70.0 73.8;;; 0.0 0.0 0.0 0.0]
    image_notreshaped=[0 0 0 0;;; 0 0 0 0;;; 0 0 0 0]
    coor =[0.5 6.24 0.0; 0.5 6.62 0.0; 0.5 7.0 0.0; 0.5 7.38 0.0]
    image=[0 0 0; 0 0 0; 0 0 0; 0 0 0]


    AAToID = Dict('M'=> Int32(1), 'N'=>Int32(3), 'A'=> Int32(2))
#        AAToID = Dict('M'=> Int32(0), 'N'=>Int32(1), 'A'=> Int32(2))

   #for (num,atom) in enumerate(set)
   #     AAToID[atom]=num
   # end


    snapshot=GSDFormat.Frame()
    snapshot.configuration.step=1
    snapshot.configuration.dimensions=3
    snapshot.configuration.box=BoxSize./10.0

    snapshot.particles.N=N
    snapshot.particles.position=coor
    IDToAA=Dict(value=>key for (key,value) in AAToID)
    snapshot.particles.types = ["M", "A", "N"]
    snapshot.particles.typeid = UInt32.([0, 2, 1, 0])

    snapshot.particles.image = image
    mass_charge=Vector{Float64}([1.0,2.0,3.0,1.0])
    snapshot.particles.mass = Vector{Float32}([1.0,2.0,3.0,1.0])
    snapshot.particles.charge = Vector{Float32}([1.0,2.0,3.0,1.0])
    snapshot.particles.diameter = [0.1, 0.3, 0.2, 0.1]

    snapshot.bonds.N=N-1
    snapshot.bonds.types = ["O-O"]
    snapshot.bonds.typeid = zeros(Int32, N-1)
    snapshot.bonds.group = [0 1; 1 2; 2 3]

    DiMap=Dict()
    DiList=Matrix{Int}(undef,4,4)

    file = GSDFormat.open(filename_test, 'w')
    GSDFormat.append(file, snapshot)
    GSDFormat.close(file)


    HPSAnalysis.Setup.writeGSDStartFile(filename,N,N-1,N-2,N-3,BoxSize,coor_notreshaped,AAToID,sequences,image_notreshaped,mass_charge,mass_charge,DiMap,DiList,AAToID,false, "Calvados2", [],[])


    frame_setup = GSDFormat.open(filename, 'r')[1]
    frame_manual = GSDFormat.open(filename_test, 'r')[1]

    for container in ["particles","bonds","angles","dihedrals","impropers","constraints","pairs"]
    #@testset "$container" begin 
        snapshot_container = getproperty(frame_setup, Symbol(container))
        test_container = getproperty(frame_manual, Symbol(container))

            for name in GSDFormat.get_container_names(test_container)
                @namedtest "$name"  getproperty(test_container,name) == getproperty(snapshot_container, name) 
                #println("$name $(getproperty(test_container,name) == getproperty(snapshot_container, name))")
            end
        #end
    end
#    @test files_are_equal(filename_test,filename)

    #rm(filename_test)
    #rm(filename)
end
