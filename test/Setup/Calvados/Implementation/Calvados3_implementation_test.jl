### Run diferent proteins for which we know the R_g from the Calvados3 paper, to verify our implementation. As multi domain proteins we use Tia1, Ubq3, Ubq2, Ubq4 and Gal3.
using Distributed
addprocs(5)

@everywhere using HPSAnalysis

PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)
ENV["PYCALL_JL_RUNTIME_PYTHON"]="$(EnvironmentPath)/bin/python"

BasePath = "$(SetupTestPath)/implementation_test/"
rm(BasePath; force=true, recursive=true)
mkpath(BasePath)

Proteins = ["tia1", "ubq2", "ubq3", "ubq4", "gal3"]
DomainDict = Dict("tia1" => [(6, 82),(95, 172),(190, 275)], "ubq2" => [(11, 82),(87, 158)], "ubq3" => [(1, 72),(77, 148),(153, 224)], "ubq4" => [(1, 72),(77, 148),(153, 224),(229, 300)], "gal3" => [(117, 250)])

ProteinJSON = Dict("tia1" => "/localscratch/test/fold_tia1/fold_tia1_full_data_0.json", "ubq2" => "/localscratch/test/fold_ubq2/fold_ubq2_full_data_0.json", "ubq3" => "/localscratch/test/fold_ubq3/fold_ubq3_full_data_0.json", "ubq4" => "/localscratch/test/fold_ubq4/fold_ubq4_full_data_0.json", "gal3" => "/localscratch/test/fold_gal3/fold_gal3_full_data_0.json")

ProteinCif  = Dict("tia1" => "/localscratch/test/fold_tia1/fold_tia1_model_0.cif", "ubq2" => "/localscratch/test/fold_ubq2/fold_ubq2_model_0.cif", "ubq3" => "/localscratch/test/fold_ubq3/fold_ubq3_model_0.cif", "ubq4" => "/localscratch/test/fold_ubq4/fold_ubq4_model_0.cif", "gal3" => "/localscratch/test/fold_gal3/fold_gal3_model_0.cif")

pH=7.0
width_multiplier=1.0
Temperatures = 300.0

@everywhere function run_sim_prot(protein, BasePath, DomainDict, ProteinJSON, ProteinCif, pH, width_multiplier, Temperatures)
    RunsPerProtein = 1
    FoldedDomains = DomainDict
    ProteinToJSON = ProteinJSON
    ProteinToCif = ProteinCif
    mkpath(BasePath*"$(protein)")
    for temp in Temperatures
        pad = RunsPerProtein
        mkpath(BasePath*"$(protein)/$(temp)K/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/Restart/")
        Path = BasePath*"$(protein)/$(temp)K/RUN_$(pad)/"
        
        println(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
        Seq = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        NChains = 1
        Sequences= [deepcopy(Seq) for _ in 1:NChains]
        Proteins = [deepcopy(protein) for _ in 1:NChains]

        Info = "SLAB Simulation script for $protein.\n\n"
        BoxLengthShort = Float32(350.0)      
        BoxLengthLong = Float32(1500.)
        BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]

        SimulName = "$(protein)_$temp"

        (pos, Data) = HPSAnalysis.CreateStartConfiguration(SimulName,Path , Float32.([BoxLengthShort,BoxLengthShort*width_multiplier , BoxLengthShort]), Proteins, Sequences, Regenerate=true; Axis="y", SimulationType="Calvados3",ProteinToDomain=FoldedDomains,ProteinToCif=ProteinToCif)

        ENM = HPSAnalysis.Setup.BuildENMModel(Data, FoldedDomains, Proteins, Sequences, ProteinToJSON)

        HPSAnalysis.Setup.writeStartConfiguration("./$(protein)_slab","./$(SimulName)_Start_slab.txt", Info, Sequences, BoxSize , 1_000, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH,domain=FoldedDomains,Device="CPU",ChargeTemperSwapSteps=10_000,WriteOutFreq=10_000, ENM)

        sim.run("$(Path)/")
    end
end

procs = workers()
i = 1
function nextproc()
    global i
    p = procs[i]
    i = i % length(procs) + 1
    return p
end
#=
@sync begin
    for protein in Proteins
        @async remotecall_wait(run_sim_prot, nextproc(), protein, BasePath, DomainDict, ProteinJSON, ProteinCif, pH, width_multiplier, Temperatures)
    end
end=#
@sync begin
    protein = "tia1"
    @async remotecall_wait(run_sim_prot, nextproc(), protein, BasePath, DomainDict, ProteinJSON, ProteinCif, pH, width_multiplier, Temperatures)
end