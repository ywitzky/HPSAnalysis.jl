


const polyply = "$(EnvironmentPath)/bin/polyply"
const martinize2 = "$(EnvironmentPath)/bin/martinize2"

module Polyply
import ..HPSAnalysis: EnvironmentPath, PkgPath

const polyply = "$(EnvironmentPath)/bin/polyply"
const martinize2 = "$(EnvironmentPath)/bin/martinize2"
using ..BioData: OneToThree
using Scratch, DataStructures,BioStructures,DSSP_jll, Printf


@doc raw"""
    ConvertGroToPDB(Path, Filename)

A gro data file is converted to a pdb data file.
    
**Arguments**
- `Path::String`: Path to the gro file.
- `Filename::String`: Name of the gro file.

**Creat**:
* A pdb data file.
"""
function ConvertGroToPDB(Path, Filename)
    GMX_Path="/localscratch/Programs/gromacs-2022.2/bin/gmx"
    run(`$GMX_Path editconf -f $(Path)$(Filename).gro -o $(Path)$(Filename).pdb`)
end

@doc raw"""
    GenerateITPFilesOfSequence(Names, Sequences, OutputPath)

Generate a itp file for one protein with the given sequence, with datas like atoms, bonds, angels.
    
**Arguments**
- `Names::Array(String)`: Names of the proteins.
- `Sequences::Array(String)`: Sequences of the proteins.
- `OutputPath::String`: Path for the itp file.

**Creat**:
* A itp file.
"""
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

@doc raw"""
    GenerateSlabTopologyFile(Filename, ITPPath, Names, SimulationName)

Generate a Topology file which include the number and name of the proteins that are simulated.
    
**Arguments**
- `Filename::String`: Path for the topology file.
- `ITPPath::String`: Pathe of the itp file.
- `Names:Array(String):`: List of protein names.
- `SimulationName::String`: Name of the simulation.

**Creat**:
* A topology data file.
"""
function GenerateSlabTopologyFile(Filename, ITPPath, Names, SimulationName)
    f = open(Filename,"w+")

    write(f,"#include \"$(PkgPath)/data/Polyply/martini_v3.0.0.itp\" \n")

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

@doc raw"""
    GenerateCoordinates(SimulationPath, SimulationName, Box, TopologyFile)

Generate a gro file from the TopologyFile which is used for the start coordinates.
    
**Arguments**
- `SimulationPath::String`: Path where to save the gro data file.
- `SimulationName::String`: Name of the simulation (save name).
- `Box::Arry(Float)`: Diameters of the simulationbox.
- `TopologyFile::String`: Topology data file.

**Creat**:
* A gro data file.
"""
function GenerateCoordinates(SimulationPath, SimulationName, Box, TopologyFile)
    max_iteration=800 #default 800
    max_force=50000.0 #default 50_000
    grid_spacing= 0.2 #default 0.2nm
    run(`$polyply gen_coords -p $TopologyFile -o $(SimulationPath)$(SimulationName).gro -name $(SimulationName) -gs $grid_spacing -mf $(max_force) -mir $(max_iteration) -box $(Box)`)
end

@doc raw"""
    readSimpleGRO(Filename, x,y,z)

Read the coordinates of the atoms from a gro file.
    
**Arguments**
- `Filename::`: Path of the pdb file.
- `x::Array(Float)`: List of x-coordinates.
- `y::Array(Float)`: List of y-coordinates.
- `z::Array(Float)`: List of z-coordinates.
"""
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

@doc raw"""
    readPDB(Filename, x,y,z)

Read the coordinates of the atoms from a pdb file.
    
**Arguments**
- `Filename::`: Path of the pdb file.
- `x::Array(Float)`: List of x-coordinates.
- `y::Array(Float)`: List of y-coordinates.
- `z::Array(Float)`: List of z-coordinates.
"""
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

@doc raw"""
    GenerateENM_ITPFilesOfSequence(BasePath::String, Names, Domains::Dict{String,Vector{Tuple{I,I}}})

Generate pdb and itp data files for the ENM with martinize2 (only Calvados3).
    
**Arguments**
- `BasePath::String`: Data path.
- `Names::List(String)`: List of protein names.
- `Domains::Dict{String,Vector{Tuple{I,I}}}`: Dictionary with domains for each protein.

**Creat**:
* Data files for the ENM.
"""
function GenerateENM_ITPFilesOfSequence(BasePath::String, Names, Domains::Dict{String,Vector{Tuple{I,I}}}) where {I<:Integer}
    Path = "$(BasePath)/InitFiles/Elastic_Files/"
    NameSet = Set(Names)
    
    for name in (NameSet)
        input_pdb_File = "$(BasePath)/InitFiles/PDBFiles/$(name).pdb"
    
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

        cd(BasePath)
        run(`$martinize2 -f $(input_pdb_File) -o $(output_pdb) -x $(cg_pdb) -ff martini3001 -ss $(ss_string) -elastic -eu 0.9 -ef $force  -eunit $domains_str -name $(name)`)

        mv("$(BasePath)/$(name)_0.itp", "$(BasePath)/InitFiles/ITPS_Files/$(name).itp"; force=true)
    end
end


end
