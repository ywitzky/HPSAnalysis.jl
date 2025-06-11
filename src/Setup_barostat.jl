function AlphaFold_startpos(ProteintoCif, Proteins, Sequences)
    length=0
    for Seq in Sequences
        length += length(Seq)
    end
    pos = zeros(Float64, length, 3)
    index = 1
    for protein in Proteins
        open(ProteintoCif[protein], "r") do io
            for line in eachline(io)
                if startswith(line, "ATOM")
                    field = split(strip(line))
                    x = parse(Float64, field[11])
                    y = parse(Float64, field[12])
                    z = parse(Float64, field[13])
                    pos[index, :] = [x, y, z] 
                end
            end
        end
    end
    return pos
end


function CreateStartConfiguration_barostat(SimulationName::String, Path::String, BoxSize::Vector{ChoosenFloatType}, Proteins::Vector{String}, Sequences::Vector{String} ; Axis="y", Regenerate=true,SimulationType="Calvados2",ProteinToDomain=Dict(),ProteinToCif=Dict())
    #Definition of Paths for the parameters
    Data = SimData()
    Data.BasePath= Path
    Data.PlotPath=Data.BasePath*"/Plots/"
    Data.DataPath=Data.BasePath*"/Data/"
    Data.xFilePath = Data.DataPath*"x.bin"
    Data.yFilePath = Data.DataPath*"y.bin"
    Data.zFilePath = Data.DataPath*"z.bin"
    Data.x_uw_FilePath = Data.DataPath*"x_uw.bin"
    Data.y_uw_FilePath = Data.DataPath*"y_uw.bin"
    Data.z_uw_FilePath = Data.DataPath*"z_uw.bin"
    Data.Reduce=1 ### only a factor for prereduced data.
    mkpath(Data.PlotPath)
    mkpath(Data.DataPath)

    #Definition of the Box
    Data.BoxLength = [Float32(BoxSize[1]),Float32(BoxSize[2]),Float32(BoxSize[3])]
    Data.BoxSize = zeros(eltype(Data.x), 3,2 )#Matrix of FloatType
    Data.BoxSize[1,1] = -Data.BoxLength[1]/2.
    Data.BoxSize[1,2] =  Data.BoxLength[1]/2.
    Data.BoxSize[2,1] = -Data.BoxLength[2]/2.
    Data.BoxSize[2,2] =  Data.BoxLength[2]/2.
    Data.BoxSize[3,1] = -Data.BoxLength[3]/2.
    Data.BoxSize[3,2] =  Data.BoxLength[3]/2.

    #Definition: Number of Steps,Chains
    Data.NSteps = 1 ### only for the creation, will be changed later.
    Data.NChains = length(Sequences)
    Data.StepFrequency=1
    Data.SimulationName=SimulationName

    Data.Sequences = Sequences

    #Definition of the length of all chains
    Data.NAtoms = 0
    for Prot in Proteins
        Data.NAtoms += length(ProteinSequences.NameToSeq[Prot]) 
    end

    #All Residues get an ID, same Residue -> same ID, and definition for other way around
    Data.IDs = zeros(eltype(Data.NAtoms), Data.NAtoms)
    cnt = Int32(1)
    ind = 1
    ResNameToID = Dict{Char, Int32}()
    for chain in Data.Sequences
        for AA in chain
            if !haskey(ResNameToID,AA)
                ResNameToID[AA] = cnt
                cnt += 1
            end
            Data.IDs[ind] = ResNameToID[AA]
            ind += 1
        end
    end
    Data.IDToResName = Dict( (v => string(k)) for (k, v) in ResNameToID)

    #Matrix of the start and end values of each chain
    Data.ChainStart = zeros(eltype(Data.NSteps), Data.NChains)
    Data.ChainStop  = zeros(eltype(Data.NSteps), Data.NChains)
    Data.ChainStart[1] = 1
    Data.ChainStop[1] = length(Data.Sequences[1])
    for (i,Seq) in enumerate(Data.Sequences)
        if i ==1 
            continue
        end
        Data.ChainStart[i] = Data.ChainStart[i-1]+length(Seq)
        Data.ChainStop[i] = Data.ChainStop[i-1]+length(Seq)
    end

    ### allocate disk space for X
    #creat path for the coordinates
    Data.xio= open(Data.xFilePath,"w+")
    Data.yio= open(Data.yFilePath,"w+")
    Data.zio= open(Data.zFilePath,"w+")
    #mmap creat 2D-Matrix for the coordinates (Lenght, Steps(=1))
    #Matrices of start coordinates (for each x,y,z-value)
    Data.x =  Mmap.mmap(Data.xio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y =  Mmap.mmap(Data.yio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z =  Mmap.mmap(Data.zio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    if Regenerate
        pos = AlphaFold_startpos(ProteinToCif, Proteins, Sequences)
        Data.x = pos[1, :]
        Data.y = pos[2, :]
        Data.z = pos[3, :]
    end

    close(Data.xio)
    close(Data.yio)
    close(Data.zio)

    ########## unfold along y-Axis so that the slab can be extended

    ### x,y,z are stored on disk and lazyly synchronised, ONLY IN READ MODE!!!
    Data.xio= open(Data.xFilePath,"r+")
    Data.yio= open(Data.yFilePath,"r+")
    Data.zio= open(Data.zFilePath,"r+")

    Data.x =  Mmap.mmap(Data.xio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y =  Mmap.mmap(Data.yio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z =  Mmap.mmap(Data.zio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    unfoldPositions(Data)

    Data.x_uw_io= open(Data.x_uw_FilePath,"r+")
    Data.y_uw_io= open(Data.y_uw_FilePath,"r+")
    Data.z_uw_io= open(Data.z_uw_FilePath,"r+")

    Data.x_uw =  Mmap.mmap(Data.x_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y_uw =  Mmap.mmap(Data.y_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z_uw =  Mmap.mmap(Data.z_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    ### periodically unwrap one Axis
    pos = zeros(eltype(Data.x), Data.NAtoms, 3)
    if Axis=="x"
        pos[:,1] .= Data.x_uw[:,1]
        pos[:,2] .= Data.y[:,1]
        pos[:,3] .= Data.z[:,1]
    elseif Axis=="y"
        pos[:,2] .= Data.y_uw[:,1]
        pos[:,1] .= Data.x[:,1]
        pos[:,3] .= Data.z[:,1]
    elseif Axis=="z"
        pos[:,3] .= Data.z_uw[:,1]
        pos[:,1] .= Data.x[:,1]
        pos[:,2] .= Data.y[:,1]
    end

    ### shift so that box center is at 0,0,0
    pos[:,1] .-= BoxSize[1]/2
    pos[:,2] .-= BoxSize[2]/2
    pos[:,3] .-= BoxSize[3]/2

    Data.x[:,1] .=  pos[:,1]
    Data.y[:,1] .=  pos[:,2]
    Data.z[:,1] .=  pos[:,3]

    Data.x_uw[:,1] .-= BoxSize[1]/2
    Data.y_uw[:,1] .-= BoxSize[2]/2
    Data.z_uw[:,1] .-= BoxSize[3]/2

    tmp_x = Data.x_uw
    tmp_z = Data.z_uw
    Data.x_uw= deepcopy(Data.x)
    Data.z_uw = deepcopy(Data.z)
    WriteAsPDB(Data; Wrapped=false)
    Data.x_uw = tmp_x
    Data.z_uw = tmp_z

    ### since they are not needed anymore
    close(Data.xio)
    close(Data.yio)
    close(Data.zio)
    close(Data.x_uw_io)
    close(Data.y_uw_io)
    close(Data.z_uw_io)

    return (pos, Data) 
end






function writeStartConfiguration_forBarostat(BasePath, fileName, StartFileName, Info, Sequences, BoxSize, ProteintoCif; NSteps=100_000_000, SimulationType="Calvados2", Temperature=300, MixingRule="1-1001-1", Pos=zeros(Float32, 0), InitStyle="AlphaFold", SaltConcentration=0.15, pH=6, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, HOOMD=false, OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma, WriteOutFreq=100_000, Device="GPU", yk_cut=40.0, ah_cut=20.0, domain=Array([]), ENM)
    
    NAtoms, NBonds, NAngles, NDihedrals, AlphaAddition, SimulationType, AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, NAtomTypes, dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list = HPSAnalysis.startConfigurationSetup(Sequences, SimulationType, pH, OneToChargeDef, OneToLambdaDef, OneToSigmaDef, MixingRule)

    if InitStyle=="AlphaFold"
        pos = AlphaFold_startpos(ProteintoCif, Proteins, Sequences)
    else
        pos = HPSAnalysis.createDenseStartingPosition(Sequences, BoxSize)
    end
    return pos
end