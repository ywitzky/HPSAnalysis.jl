
module Polyply
using HPSAnalysis
using ..BioData: OneToThree
using Scratch, DataStructures,BioStructures,DSSP_jll, Printf


polyply = "$(HPSAnalysis.EnvironmentPath)/bin/polyply"
martinize2 = "$(HPSAnalysis.EnvironmentPath)/bin/martinize2"


polyply="/localscratch/Python3.13Environments/OldPython/bin/polyply"

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
       run(`$polyply gen_params -name $(name) -lib martini3 -seq $SeqString -o $(OutputPath)$(name).itp`)
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
        write(f, "$(name)_0 $(Occurences[name])\n")
    end
    close(f)
end

function GenerateCoordinates(SimulationPath, SimulationName, Box, TopologyFile)
    max_iteration=800 #default 800
    max_force=50000.0 #default 50_000
    grid_spacing= 0.2 #default 0.2nm
    run(`$polyply gen_coords -p $TopologyFile -o $(SimulationPath)$(SimulationName).gro -name $(SimulationName) -gs $grid_spacing -mf $(max_force) -mir $(max_iteration) -box $(Box)`)
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

### Generate ITP files with AlphaFold Based Elastic Network model - not the one taken for the final sim
function GenerateENM_ITPFilesOfSequence(Sim::HPSAnalysis.SimData{T,I},Names,  Domains::Dict{String,Vector{Tuple{I2,I2}}}) where {T<:Real, I<:Integer, I2<:Integer}
    Path = "$(Sim.BasePath)/InitFiles/Elastic_Files/"
    NameSet = Set(Names)
    
    for name in (NameSet)
        input_pdb_File = "$(Sim.BasePath)/InitFiles/PDBFiles/$(name).pdb"
    
        ### compute secondary structure assigment according to DSSP
        struc = read(input_pdb_File, BioStructures.PDBFormat, run_dssp=true)
        file=sscode.(collectresidues(struc))

        ss_string_ = join([c == ' ' ? 'C' : c for c in file])
        ### turn P and missing into C (since martinize cant parse this info)
        ss_string = replace(join(file), 'P' => 'C', ' ' => 'C')

        domains_str=join(map(r->"$(r[1]):$(r[2])", Domains[name]), ",")

        output_pdb = "$(Path)$(name).top"
        cg_pdb = "$(Path)$(name)_cg_protein.pdb"

        force=4000 ## default 700; 4000 value of back bone bonds


        run(`$martinize2 -f $(input_pdb_File) -o $(output_pdb) -x $(cg_pdb) -ff martini3001 -ss $(ss_string) -b /uni-mainz.de/homes/ywitzky/.julia/scratchspaces/f54fdea9-a25b-4801-be39-89e20f64ecdd/test/Setup/RS31/300K/RUN_001/InitFiles/build_file.bld -elastic -eu 0.9 -ef $force  -eunit $domains_str -name $(name)`)

        mv("$(Sim.BasePath)/$(name)_0.itp", "$(Sim.BasePath)/InitFiles/ITPS_Files/$(name).itp"; force=true)
    end
end

### Rewrite AlphaFold Cif InputFiles to PDB in Folder
function RewriteCifToPDB(Sim::HPSAnalysis.SimData{T,I}, ProteinToCif,Proteins) where {T<:Real, I<:Integer}
    subpath = "$(Sim.BasePath)/InitFiles/PDBFiles"
    cifpath = "$(Sim.BasePath)/InitFiles/CifFiles"
    mkpath(subpath)
    mkpath(cifpath)
    for Prot in Set(Proteins)
        CifPath = ProteinToCif[Prot]
        try
            readlines(CifPath)
        catch
            @error "There is no AlphaFold data for protein $(Prot) at $(CifPath)."
        end

        cp(CifPath, "$(cifpath)/$(Prot).cif"; force=true)
        pdb_file = open("$(subpath)/$(Prot).pdb", "w")
        for line in readlines(CifPath)
            fields = strip.(split(line))
            if !isempty(fields) && fields[1] == "ATOM"
                ID = parse(Int, fields[2])
                symbole = fields[3]
                label_atom = fields[4]
                comp = fields[6]
                asym = fields[7]
                seq  = parse(Int, fields[9])
                x = parse(Float64, fields[11])
                y = parse(Float64, fields[12])
                z = parse(Float64, fields[13])
                iso_or_equiv= parse(Float64, fields[15])

                println(pdb_file, @sprintf("ATOM  %5d %4s %-3s %s%4d    %8.3f%8.3f%8.3f  1.00  %1.2f %10s", ID, label_atom, comp, asym, seq, x, y, z, iso_or_equiv, symbole))
            end
        end
        close(pdb_file)
    end
end

end
