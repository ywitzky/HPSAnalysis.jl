
module Polyply
using HPSAnalysis
using ..BioData: OneToThree
using Scratch, DataStructures
polyply = "$(HPSAnalysis.EnvironmentPath)/bin/polyply"


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
        run(`$polyply gen_params  -name $(name) -seq $SeqString -lib martini3  -o $(OutputPath)$(name).itp`) 
    end
end

function GenerateSlabTopologyFile(Filename, ITPPath, Names, SimulationName)
    f = open(Filename,"w+")

    write(f,"#include \"$(HPSAnalysis.PkgPath)/data/Polyply/martini_v3.0.0.itp\" \n")

    write(f, "[ atomtypes ]\n")
    write(f, "VS 0.00 0.000 V 0.0 0.0\n\n")

    write(f, "[ nonbond_params ]\n")
    write(f, "VS    W   1 0.4650000000    0.5000000000\n\n")
    
    NameSet = Set(Names)
    Occurences = counter(Names)
    for name in NameSet
        write(f, "#include \"$(ITPPath)$(name).itp\" \n")
    end
    write(f,"\n[ system ]\n; name\nProtein Sim\n\n[ molecules ]\n; name  number\n")

    for (ind, name) in enumerate(NameSet)
        write(f, "$(name) $(Occurences[name])\n")
    end
    close(f)
end

function GenerateCoordinates(SimulationPath, SimulationName, Box, TopologyFile)
    run(`$polyply gen_coords -p $TopologyFile -o $(SimulationPath)$(SimulationName).gro -name $(SimulationName) -box $(Box)`)
end

function readSimpleGRO(Filename, x,y,z)
    f = open(Filename)
    readline(f) ### command that generate the file
    readline(f) ### the number of lines
    ind=1
    while !eof(f) 
        line = readline(f)
        if length(line)>15 && line[14:15]=="BB"#"CA"
            x[ind, 1] = parse(Float32, line[21:28])
            y[ind, 1] = parse(Float32, line[29:37])
            z[ind, 1] = parse(Float32, line[38:44])
            ind+=1
        end
    end
    ### convert from nm to angstroem
    x .*= 10.0
    y .*= 10.0
    z .*= 10.0

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
#end

using BioStructures
using DSSP_jll
function ElasticNetworkModel(Names, Sequences, Path::String,pdb_Path, domains)
    Path = "$(Path)/Elastic_Files/"
    NameSet = Set(Names)
    #println(NameSet)
    martinize2 = "$(HPSAnalysis.EnvironmentPath)/bin/martinize2"
    input_pdb_File = pdb_Path
    #ss_file = "$(Path)ss.txt"
    
    struc = read(input_pdb_File, BioStructures.PDBFormat, run_dssp=true)
    file=sscode.(collectresidues(struc))
    ss_string = join([c == ' ' ? 'C' : c for c in file])
    ### turn P into C (because of BioStructures and martinize)
    ss_string=replace(ss_string, 'P' => 'C')

    #println(ss_string)
    #write(ss_file, ss_string)

    #[[1,70],[90,150]]=>"1:70,90:150"
    domains_str=join(map(r->"$(r[1]):$(r[2])", domains), ",")
    #print(domains_str)
    for name in (NameSet)
        output_pdb = "$(Path)$(name).top"
        cg_pdb = "$(Path)$(name)_cg_protein.pdb"

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

        #run(`$polyply gen_params  -name $(name) -seq $SeqString -lib martini3  -o $(OutputPath)$(name).itp`)
        run(`$martinize2 -f $(input_pdb_File) -o $(output_pdb) -x $(cg_pdb) -ff martini3001 -ss $(ss_string) -elastic -eunit $domains_str -name $(name)`)
    end
    if false
        #martinize2 = "$(HPSAnalysis.EnvironmentPath)/bin/martinize2"

        # File paths
        #cifFile = "/localscratch/lafroehl/Hiwi/fold_rs31/test.pdb"
        #output_pdb = "$(Path).top"
        #cg_pdb = "$(Path)_cg_protein.pdb"
        #ss_file = "$(Path)_ss.txt"

        struc = read(cifFile, BioStructures.PDBFormat, run_dssp=true)
        file=sscode.(collectresidues(struc))
        ss_string = join([c == ' ' ? 'C' : c for c in file])
        ### turn P into C (because of BioStructures and martinize)
        ss_string=replace(ss_string, 'P' => 'C')

        #println(ss_string)
        write(ss_file, ss_string)

        #[[1,70],[90,150]]=>"1:70,90:150"
        domains_str=join(map(r->"$(r[1]):$(r[2])", domains), ",")
        #print(domains_str)
        
        #creat a file with position and connection information
        run(`$martinize2 -f /localscratch/lafroehl/Hiwi/fold_rs31/test.pdb -o $(Path)_test.top -x $(Path)_cg_protein.pdb -ff martini3001 -ss $(ss_string) -elastic -eunit $domains_str`)
    end
end

using Printf
function folded_data(Proteins, Sequences, SimPath, Path, to_pdb,domains)
    try
        readlines(Path)
    catch
        @warn "The are no AlphaFold data at $(Path)"
        return
    end

    AlphaFold = readlines(Path)
    pdb_lines = String[]

    for line in AlphaFold
        fields = strip.(split(line))
        if !isempty(fields) && fields[1] == "ATOM"
            ID = parse(Int, fields[2])
            symbole = fields[3]
            label_atom = fields[4]
            #if label_atom in ["N","CA","C","O"]
            #    label_atom = "BB"
            #else
            #    label_atom = "S1"
            #end
            comp = fields[6]
            asym = fields[7]
            seq  = parse(Int, fields[9])
            x = parse(Float64, fields[11])
            y = parse(Float64, fields[12])
            z = parse(Float64, fields[13])
            iso_or_equiv= parse(Float64, fields[15])

            pdb_line = @sprintf("ATOM  %5d %4s %-3s %s%4d    %8.3f%8.3f%8.3f  1.00  %1.2f %10s", ID, label_atom, comp, asym, seq, x, y, z, iso_or_equiv, symbole)
            push!(pdb_lines, pdb_line)
        end
    end
    open(to_pdb, "w") do io
        for line in pdb_lines
            println(io, line)
        end
    end
    ElasticNetworkModel(Proteins, Sequences, SimPath, to_pdb,domains)
end

end