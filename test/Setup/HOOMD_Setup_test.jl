using HPSAnalysis,Test,Printf

filename_test="$SetupTestPath/HOOMD_write_test.csv"
filename="$SetupTestPath/HOOMD_write.csv"

@testset "WriteHOOMDSequences" begin
    sequences=["MRPVFV","MRPVF","MRPV","MRP"]
    io=open(filename_test,"w")
    for (_,seq) in enumerate(sequences)
        write(io,"$seq\n")
    end
    close(io)
    
    HPSAnalysis.Setup.WriteHOOMDSequences(filename,sequences)
    @test files_are_equal(filename_test,filename)
end

@testset "WriteDictionaries" begin
    io=open(filename_test,"w")
    ToCharge_test=Dict('A'=>1,'B'=>-1,'C'=>2,'D'=>-2)
    ToID_test=Dict('A'=>1,'B'=>2,'C'=>3,'D'=>4)
    ToMass_test=Dict('A'=>1,'B'=>2,'C'=>3,'D'=>4)
    ToDiameter_test=Dict('A'=>1,'B'=>2,'C'=>3,'D'=>4)
    ToLambda_test=Dict('A'=>1,'B'=>2,'C'=>3,'D'=>4)

    write(io, "// ID,resname, Charge, Mass, λ   \n")
    keys=['A','B','C','D']
    for key in keys
        write(io, " $(ToID_test[key]), $key, $(ToCharge_test[key]), $(ToMass_test[key]), $(ToDiameter_test[key]), $(ToLambda_test[key])\n")
    end

    close(io)
    HPSAnalysis.Setup.WriteDictionaries(filename,ToCharge_test,ToID_test,ToMass_test,ToDiameter_test,ToLambda_test)
    @test files_are_equal(filename_test,filename)
end

@testset "WriteHOOMDParticlesInput" begin
    sequences=["MRPVFV","MRPVF","MRPV","MRP"]
    io=open(filename_test,"w")
    set=Set(join(sequences))
    To=Dict()
    coor=fill(0.5,4,maximum(length(seq) for seq in sequences),3)
    for (num,atom) in enumerate(set)
        To[atom]=num
    end
    ToCharge=To
    ToID=To
    ToMass=To
    ToDiameter=To
    pos=coor
    image=coor

    write(io, "### N, id, x ,y ,z , charge, m, diameter \n");
    cnt=1
    for (SeqID, Seq) in enumerate(sequences)
        for (atom,res) in enumerate(Seq)
            write(io, "$(cnt), $(ToID[res] -1) , $(@sprintf("%.3f",pos[SeqID,atom, 1])), $(@sprintf("%.3f",pos[SeqID,atom, 2])), $(@sprintf("%.3f",pos[SeqID,atom, 3])), $(ToCharge[res]), $(ToMass[res]), $(ToDiameter[res]) , $(image[SeqID, atom,1]) , $(image[SeqID, atom,2]) , $(image[SeqID, atom,3]) \n");
            cnt += 1
        end
    end

    close(io)
    HPSAnalysis.Setup.WriteHOOMDParticlesInput(filename, pos, ToCharge, ToID, sequences, ToMass, ToDiameter, image)
    @test files_are_equal(filename_test,filename)
end


@testset "WriteParams" begin
    sequences=["MRPVFV","MRPVF","MRPV","MRP"]
    io=open(filename_test,"w")
    SimName="C2"; Temp=300; NSteps=10000; NOut=10; Timestep=1000
    Box=Array([10,101,10])
    Seed=10; Minimise=true; TrajectoryName="testtraj.gsd"; UseAngles=true; UseCharge=true; Alt_GSD_Start="-"; Create_Start_Config=false; ϵ_r=1.73136; κ=1.0; Device="GPU"; yk_cut=4.0; ah_cut=2.0; ionic=0.1; pH=7.0;SimType="Calvados2";domain=Array([[0,0]]);
    
    write(io, "Simname: $SimName\n")
    write(io, "Domains: $(domain)\n")
    write(io, "Seed: $Seed\n")
    write(io, "Temp: $Temp\n")
    write(io, "ionic: $ionic\n")
    write(io, "pH: $pH\n")
    write(io, "SimulationType: $SimType\n")
    write(io, "NSteps: $NSteps\n")
    write(io, "NOut: $NOut\n")
    write(io, "dt: $Timestep\n")
    write(io, "Lx: $(Box[1])\n")
    write(io, "Ly: $(Box[2])\n")
    write(io, "Lz: $(Box[3])\n")
    write(io, "Minimise: $(Minimise)\n")
    write(io, "Trajectory: $(TrajectoryName)\n")
    write(io, "UseAngles: $(UseAngles)\n")
    write(io, "UseCharge: $(UseCharge)\n")
    write(io, "Alt_GSD_Start: $(Alt_GSD_Start)\n")
    write(io, "Create_Start_Config: $(Create_Start_Config)\n")
    write(io, "epsilon_r: $(ϵ_r)\n")
    write(io, "yk_prefactor: $(138.9315360433804/ϵ_r)\n")
    write(io, "kappa: $(κ)\n")
    write(io, "Device: $(Device)\n")
    write(io, "YukawaCutoff: $(yk_cut)\n")
    write(io, "AHCutoff: $(ah_cut)\n")
    close(io)
    
    HPSAnalysis.Setup.WriteParams(filename, SimName, Temp, NSteps, NOut, Timestep, Box, Seed; Minimise, TrajectoryName, UseAngles, UseCharge, Alt_GSD_Start, Create_Start_Config, ϵ_r, κ, Device, yk_cut, ah_cut, ionic, pH)
    @test files_are_equal(filename_test,filename)
end
