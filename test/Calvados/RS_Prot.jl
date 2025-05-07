using HPSAnalysis

PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)
ENV["PYCALL_JL_RUNTIME_PYTHON"]="$(EnvironmentPath)/bin/python"

using PyCall

#rm(SetupTestPath; force=true, recursive=true)
mkpath("$SetupTestPath/HOOMD_Setup/")
BasePath=SetupTestPath


ToCreate =  ["RS31"]
FoldedDomains = Dict("RS31" => [(1,70),(90,155)])
ProteinToCif =Dict("RS31" => "/localscratch/test/fold_rs31/fold_RS31_model_0.cif")
ProteinToJSON =Dict("RS31" => "/localscratch/test/fold_rs31/fold_rs31_full_data_0.json")


Temperatures=300
RunsPerProtein=1
pH=7.0
#SideLength=0.75
width_multiplier=1.0

pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
sim = pyimport("Submit_HOOMD")

for (protID, protein) in enumerate(ToCreate)
    mkpath(BasePath*"$(protein)/")
    for temp in Temperatures
        pad = "001"
        mkpath(BasePath*"$(protein)/$(temp)K/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/Restart/")
        Path = BasePath*"$(protein)/$(temp)K/RUN_$(pad)/"
        cd(Path)

        Seq = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        NChains=4

        Seq = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        Sequences= [deepcopy(Seq) for _ in 1:NChains]
        Proteins = [deepcopy(protein) for _ in 1:NChains]

        ###FoldedDomain -> NChains * FoldedDomain
        #=
        NewDomain = copy(FoldedDomain)
        for i in 1:NChains-1
            shift = i * length(Seq)
            for dom in FoldedDomain
                push!(NewDomain, [dom[1]+shift, dom[2]+shift])
            end
        end=#

        Info ="SLAB Simulation script for $protein.\n\n"
        BoxLengthShort=Float32(350.0)
        BoxLengthLong=Float32(1500.)
        BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]

        SimulName = "$(protein)_$temp"

        (pos, Data) = HPSAnalysis.CreateStartConfiguration(SimulName,Path , Float32.([BoxLengthShort,BoxLengthShort*width_multiplier , BoxLengthShort]), Proteins, Sequences, Regenerate=false; Axis="y", SimulationType="Calvados3",ProteinToDomain=FoldedDomains,ProteinToCif=ProteinToCif)

        ENM = HPSAnalysis.Setup.BuildENMModel(Data, FoldedDomains, Proteins, Sequences, ProteinToJSON)

        itp_Path = "$(Data.BasePath)/InitFiles/ITPS_Files/$(protein).itp"
        #=
        HPSAnalysis.Setup.writeStartConfiguration("./$(protein)_slab","./$(SimulName)_Start_slab.txt", Info, Sequences, BoxSize , 30_000, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH,domain=NewDomain,Device="CPU",ChargeTemperSwapSteps=10_000,WriteOutFreq=10_000, itp_Path)

        sim.run("$(Path)/")
        =#
    end
end
