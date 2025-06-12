using HPSAnalysis
using PyCall

BasePath = "$(SetupTestPath)/Barostat_test/"
if isdir(BasePath)
    rm(BasePath; force=true, recursive=true)
end
mkpath(BasePath)

Proteins = ["tia1", "ubq2", "ubq3", "ubq4", "gal3"]
DomainDict = Dict("tia1" => [(6, 82),(95, 172),(190, 275)], "ubq2" => [(11, 82),(87, 158)], "ubq3" => [(1, 72),(77, 148),(153, 224)], "ubq4" => [(1, 72),(77, 148),(153, 224),(229, 300)], "gal3" => [(117, 250)])

ProteinJSON = Dict("tia1" => "/localscratch/test/fold_tia1/fold_tia1_full_data_0.json", "ubq2" => "/localscratch/test/fold_ubq2/fold_ubq2_full_data_0.json", "ubq3" => "/localscratch/test/fold_ubq3/fold_ubq3_full_data_0.json", "ubq4" => "/localscratch/test/fold_ubq4/fold_ubq4_full_data_0.json", "gal3" => "/localscratch/test/fold_gal3/fold_gal3_full_data_0.json")

ProteinCif  = Dict("tia1" => "/localscratch/test/fold_tia1/fold_tia1_model_0.cif", "ubq2" => "/localscratch/test/fold_ubq2/fold_ubq2_model_0.cif", "ubq3" => "/localscratch/test/fold_ubq3/fold_ubq3_model_0.cif", "ubq4" => "/localscratch/test/fold_ubq4/fold_ubq4_model_0.cif", "gal3" => "/localscratch/test/fold_gal3/fold_gal3_model_0.cif")

pH = 7.0
width_multiplier = 1.0
Temperatures = 300.0
function run_sim_prot(protein, BasePath, DomainDict, ProteinJSON, ProteinCif, pH, width_multiplier, Temperatures)
    PkgSourcePath = "/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
    EnvironmentPath = HPSAnalysis.getPythonEnvironment(PkgSourcePath)
    ENV["PYCALL_JL_RUNTIME_PYTHON"] = "$(EnvironmentPath)/bin/python"
    pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
    sim = pyimport("Submit_HOOMD")
    func = pyimport("PythonFuncs")

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
        
        Seq = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        NChains = 6
        Sequences = [deepcopy(Seq) for _ in 1:NChains]
        Proteins  = [deepcopy(protein) for _ in 1:NChains]

        Info = "SLAB Simulation script for $protein.\n\n"
        BoxLengthShort = Float32(350.0)      
        BoxLengthLong  = Float32(1500.)
        BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]

        SimulName = "$(protein)_$temp"

        (pos, Data) = HPSAnalysis.CreateStartConfiguration_barostat(SimulName,Path , Float32.([BoxLengthShort,BoxLengthShort*width_multiplier , BoxLengthShort]), Proteins, Sequences, Regenerate=true; Axis="y", SimulationType="Calvados3",ProteinToDomain=FoldedDomains,ProteinToCif=ProteinCif)
        HPSAnalysis.RewriteCifToPDB(Path, ProteinToCif, Proteins)

        #pos = HPSAnalysis.writeStartConfiguration_Barostat(BasePath, fileName="$protein", StartFileName="$protein", Info, Sequences, BoxSize, ProteinCif)
        ENM = HPSAnalysis.Setup.BuildENMModel(Data, FoldedDomains, Proteins, Sequences, ProteinToJSON)

        HPSAnalysis.Setup.writeStartConfiguration(Path, "/$(protein)_slab","/$(SimulName)_Start_slab", Info, Sequences, BoxSize , 10_000, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH,domain=FoldedDomains,Device="CPU",WriteOutFreq=1_000, ENM)

        #sim.run("$(Path)/")
    end
end

run_sim_prot(Proteins[1], BasePath, DomainDict, ProteinJSON, ProteinCif, 7.0,0.8,300)