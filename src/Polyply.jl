
module Polyply
#include("BioData.jl")
#include("ProteinSequences.jl")
#using .BioData, .ProteinSequences
using ..BioData: OneToThree
using DataStructures

function ConvertGroToPDB(Path, Filename)
    GMX_Path="/localscratch/Programs/gromacs-2022.2/bin/gmx"
    run(`$GMX_Path editconf -f $(Path)$(Filename).gro -o $(Path)$(Filename).pdb`)
end

function GenerateITPFilesOfSequence(Names, Sequences, OutputPath)
    NameSet = Set(Names)
    for name in (NameSet)
        ind = findfirst(x-> x==name, Names)
        SeqString = []
        for AA in Sequences[ind]
            ### treat phosphorylated resdidues like 
            if AA == '#'
                push!(SeqString, "$(OneToThree['S']):1")
            elseif AA=='&'
                push!(SeqString, "$(OneToThree['T']):1")
            elseif AA=='@'
                push!(SeqString, "$(OneToThree['Y']):1")
            else
                push!(SeqString, "$(OneToThree[AA]):1")
            end
        end
       run(`polyply gen_params -name $(name) -lib martini3 -seq $SeqString -o $(OutputPath)$(name).itp`)
    end
end

function GenerateSlabTopologyFile(Filename, ITPPath, Names, SimulationName)
    f = open(Filename,"w+")

    write(f,"#include \"/localscratch/Programs/Polyply/martini_v300/martini_v3.0.0.itp\" \n")
    NameSet = Set(Names)
    Occurences = counter(Names)
    for name in NameSet
        write(f, "#include \"$(ITPPath)$(name).itp\" \n")
    end
    write(f,"\n[ system ]\n; name\nProtein Sim\n\n[ molecules ]\n; name  number\n")

    for (ind, name) in enumerate(NameSet)
        write(f, "$name $(Occurences[name])\n")
    end
    close(f)
end

function GenerateCoordinates(SimulationPath, SimulationName, Box, TopologyFile)
    run(`polyply gen_coords -p $TopologyFile -o $(SimulationPath)$(SimulationName).gro -name $(SimulationName) -box $Box`)
end

function readPDB(Filename, x,y,z)
    f = open(Filename)
    ind =1
    while !eof(f) 
        line = readline(f)
        if length(line)<6
            continue
        elseif line[1:6]=="ATOM  " 
            if line[14:15]=="BB"
                x[ind, 1] = parse(Float32, line[31:38])
                y[ind, 1] = parse(Float32, line[39:46])
                z[ind, 1] = parse(Float32, line[47:54])
                ind += 1
            end
        end
    end
end

end #Polyply

#=

SimulationName = "TestSimulation"
BoxSize=[15, 15,70] 
WorkPath="/localscratch/Programs/Polyply/"
Proteins=["RS31","RS31","RS31","RS31","RS31a","RS40", "RS41" ]
Proteins = []
Prot = "RS41"
for i in 1:110
    push!(Proteins, Prot)
end
Seqs = [] #eros(String, length(Proteins))
for (ind, Prot) in enumerate(Proteins)
    push!(Seqs , ProteinSequences.NameToSeq[Prot])
end

### generate Martini ITP Files
#GenerateITPFilesOfSequence(Proteins, Seqs, "$(WorkPath)ITPS_Files/")

### Generate Topology files
TopologyFile = "$(WorkPath)TestTopology.top"
#GenerateSlabTopologyFile(TopologyFile,"$(WorkPath)ITPS_Files/", Proteins, SimulationName)

### generate coordinates
#GenerateCoordinates(WorkPath, SimulationName, BoxSize, TopologyFile)

#convert to PDB
#ConvertGroToPDB(WorkPath, SimulationName)


NumAtoms = zeros(Int32, 1)
for Prot in Proteins
    NumAtoms[1] = length(ProteinSequences.NameToSeq[Prot]) + NumAtoms[1]
end
NAtoms = NumAtoms[1]

x  = zeros( (NAtoms, 1))
y  = zeros( (NAtoms, 1))
z  = zeros( (NAtoms, 1))

readPDB("$(WorkPath)$SimulationName.pdb", x,y,z)
=#