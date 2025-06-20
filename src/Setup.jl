include("./Setup/Polyply.jl")

module Setup

export createDenseStartingPosition, writeCollectedSlurmScript

using GSDFormat, HPSAnalysis, JSON, Printf

include("../data/BioData.jl")
include("Setup/HOOMD_Setup.jl")

using .BioData


function determineDihedrals(Sequences, Types, TypeToId, OneToHPSDihedral0110, OneToHPSDihedral1001, MixingRule="1-1001-1")
    dihedral_short_map = Dict()
    dihedral_long_map = Dict()

    dihedral_eps = ones(400)*Inf
    dihedral_cnt=0
    eps=0.
    AA_ind = 1
    AA2_ind = 1

    ReverseDihedralMap = Dict()
    key=()
    M = MixingRule=="1-1001-1" ? 4 : 2
    dihedral_list = zeros(Int32, sum(length.(Sequences).-3), M)
    index = 1
    for Sequence in Sequences
        SeqLength = length(Sequence)
        for (ind, AA) in enumerate(Sequence)
            if ind> (SeqLength-3) 
                continue
            else
                if MixingRule=="1001"
                    AA2 = Sequence[ind+3]
                    AA_ind = TypeToId[AA]
                    AA2_ind = TypeToId[AA2]
                    eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2])/2.
                    key = (min(AA_ind,AA2_ind), max(AA_ind,AA2_ind))
                elseif MixingRule=="1-1001-1"
                    if (ind ==1 )
                        AAminus = "Z"
                        AA_min_ind = 0
                    else
                        AAminus = Sequence[ind-1]
                        AA_min_ind = TypeToId[AAminus]
                    end
                    if ( ind <=(SeqLength-4) ) 
                        AA3= Sequence[ind+4]
                        AA3_ind = TypeToId[AA3]
                    else
                        AA3="X"
                        AA3_ind =-1
                    end
                    AA2 = Sequence[ind+3]
                    AA_ind = TypeToId[AA]
                    AA2_ind = TypeToId[AA2]
                    if (ind ==1)
                        eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2]+OneToHPSDihedral1001[AA3])/3.
                    elseif ( ind >(SeqLength-4) ) 
                        eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2]+OneToHPSDihedral1001[AAminus])/3.
                    else
                        eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2]+OneToHPSDihedral1001[AAminus]+OneToHPSDihedral1001[AA3])/4.
                    end
                    key = (sort([AA_min_ind,AA_ind, AA2_ind, AA3_ind]))
                elseif MixingRule=="0110"
                    AA1 = Sequence[ind+1]
                    AA2 = Sequence[ind+2]
                    AA_ind = TypeToId[AA1]
                    AA2_ind = TypeToId[AA2]
                    eps = (OneToHPSDihedral0110[AA1]+OneToHPSDihedral0110[AA2])/2.
                    key = (min(AA_ind,AA2_ind), max(AA_ind,AA2_ind))
                end
                if !haskey(dihedral_short_map, key)
                    eps_short = round(eps;digits=4)
                    same_eps_ind = findfirst(x->x==eps_short, dihedral_eps)
                    if same_eps_ind===nothing
                        dihedral_cnt+=1
                        dihedral_short_map[key]=dihedral_cnt
                        dihedral_long_map[key]=dihedral_cnt

                        ReverseDihedralMap[eps_short] = key
                        dihedral_eps[dihedral_cnt] = eps_short
                    else
                        dihedral_long_map[key]=same_eps_ind
                        key = ReverseDihedralMap[eps_short]
                    end
                end
                dihedral_list[index, :] .= key 
                index += 1
            end
        end
    end
    dihedral_eps = dihedral_eps[1:dihedral_cnt]
    return (dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list)
end

BondStrength = Dict("Calvados2"=>19.19,"ArashModell"=>20.0 )
EpsilonAshbaughHatch = Dict("Calvados2"=> 0.20,"ArashModell"=>0.2)
Cutoff = Dict("Calvados2"=> 40.0, "HPS-Alpha"=>35.0,"ArashModell"=>35.0 )

Formats=Dict([("h5md", "position image create_group yes")])

function writeHPSLammpsScript(fileName, StartFileName, AtomTypes, LongAtomTypes, AaToId, LongAtomTypesToRes, OneToCharge, OneToSigma, OneToLambda, Dihedral_eps, InitStyle, SimulationType, Temperature, AlphaStruct, Restart, NTimeSteps; WriteOutFreq=100_000,SaltConcentration=-1, pH = 6, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, OutFormat="h5md")

    ChargeTemperSim=length(ChargeTemperSteps)>0

    (ϵ_r, κ) = DetermineYukawaInteractions(;SimulationType=SimulationType, Temperature=Temperature, SaltConcentration=SaltConcentration)

    file = open(fileName, "w")
    write(file,"
    ### Initialisation

    units real ### https://docs.lammps.org/stable/units.html
    dimension 3
    newton on ### might be faster if off depends on parallelization https://docs.lammps.org/stable/newton.html
    processors * * * grid numa #map zyx
    ### processors pre dimension, can optimize a lot here.
    boundary p p p # periodic in each dimension 
    atom_style full 

    atom_modify id yes  #ids will be assigned to each atom 
    #atom_modify map hash 
    #atom_modify 10000 2 #resort ids every 10000 steps for cache improvements 

    # Reading file i easier\n")


    if ChargeTemperSim
        write(file,"\n    ### set up H-REMD
    variable q world $(join(string.(ChargeTemperSteps).*" "))
    variable id world $(join(string.(collect(1:length(ChargeTemperSteps))).*" "))
    variable energyfile string \"energy_\${id}.dat\" \n")
    else 
        write(file," variable energyfile string \"energy.dat\" \n")
    end

    if (Restart)
    write(file, "read_restart ./sim$(ChargeTemperSim ?  "_\${id}" : "").restart \n")
    else
    write(file, "read_data $(StartFileName)$(ChargeTemperSim ?  "_\${id}" : "").txt\n")
    end

    write(file, "neighbor 3.5 multi\n")
    if SimulationType=="HPS-Alpha"
        write(file,"neigh_modify every 5 delay 0
    neigh_modify collection/interval 5 20 22 24 26 36\n")
    elseif SimulationType=="Calvados2"
        write(file,"neigh_modify every 5 delay 0
    neigh_modify collection/interval 2 20.5 40.5\n")
    end

    write(file,"
    bond_style  harmonic")
    if AlphaStruct
        write(file,"
    angle_style bch
    dihedral_style gaussian")
    end
    write(file,"

    dielectric  $(@sprintf("%.5f", ϵ_r))
    
    # bond potential parameters
    bond_coeff          1   $(@sprintf("%.5f", get!(BondStrength, SimulationType,10.000)))    3.800000") 

    if AlphaStruct
        write(file,"

        # angle potential parameters 
        angle_coeff         1    4.300000

        # dihedral angle parameters\n")
        for (i, eps) in enumerate(Dihedral_eps)
            write(file,"dihedral_coeff        $(i)    $(eps)\n")
        end
    end
    write(file,"\n
    ### Coulomb and Ashbaugh Hatch Potential
    pair_style  ljlambda $(@sprintf("%.5f", 1/κ)) 0.0 $(Cutoff[SimulationType])
    dielectric  $(@sprintf("%.5f", ϵ_r))\n")
    if AlphaStruct
        write(file,"special_bonds lj/coul 0.0 0.0 0.0\n")
    else
        write(file,"special_bonds lj/coul 0.0 1.0 1.0\n")
    end
    #    special_bonds ljlambda 0.0 1.0 1.0\n")

    #for (id, ResType) in enumerate(AtomTypes)
    #    write(file, "pair_coeff $id $id $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType])  $(BioData.OneToHPSCharge[ResType])\n")
    #end



    for  ResType in LongAtomTypes
        id = AaToId[ResType]
        for  ResType2 in LongAtomTypes
            id2 = AaToId[ResType2]

            if id2<id
                continue
            end
            charge = OneToCharge[ResType]*OneToCharge[ResType2]
            lambda = (OneToLambda[ResType]+OneToLambda[ResType2])/2.
            sigma = (OneToSigma[ResType]+OneToSigma[ResType2])/2.
            if SimulationType=="HPS-Alpha" ||  SimulationType=="ArashModell"
                cutoff=4*sigma
                if abs(charge) >0.001 
                    cutoff2= Cutoff[SimulationType]
                else
                    cutoff2= 0
                end
            elseif SimulationType=="Calvados2"
                cutoff=20
                if abs(charge) >10^-8
                    cutoff2= 40
                else
                    cutoff2= 0
                end
            end
            write(file, "pair_coeff $id $id2  $(@sprintf("%.5f", get!(EpsilonAshbaughHatch, SimulationType,0.2000))) $(@sprintf("%.5f", sigma)) $(@sprintf("%.5f", lambda)) $(cutoff) $(cutoff2)\n")
        end
        #write(file, "pair_coeff $id $id ah $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType]) \n")
    end
    #write(file, "\npair_modify mix arithmetic shift yes\n")
    Charged = ""
    Uncharged = ""
    for (id, ResType) in enumerate(LongAtomTypes)
        if OneToCharge[ResType]!=0
            Charged *= "$id "
        else
            Uncharged *= "$id "
        end
    end
    if length(Charged)>0 && length(Uncharged)>0
    write(file, "group charged type $Charged
        group uncharged type $Uncharged
        
        comm_style tiled
        fix 4 all balance 100000 1.1 rcb weight group 2 charged 8 uncharged 1
        balance 1.1 rcb weight group 2 charged 8 uncharged 1\n        ")
    end

    if InitStyle=="Slab"
        
        write(file, "#comm_style tiled
    #fix 4 all balance 50000 1.05 rcb\n")
    end

    
    if ~(Restart)
        write(file,"\nminimize 1.0e-4 1.0e-6 500 1000\n")
        write(file, "reset_timestep 0\n")
    end

    write(file, "run_style verlet\n")

    if OutFormat=="xyz"
        write(file, "dump trajectory all xyz $WriteOutFreq ./Trajectory$(ChargeTemperSim ?  "_\${id}" : "").xyz")
    elseif OutFormat=="xtc"
        write(file,"dump trajectory all xtc $WriteOutFreq ./Trajectory$(ChargeTemperSim ?  "_\${id}" : "").xtc\n")
    elseif OutFormat=="h5md"
        write(file,"dump trajectory all h5md $WriteOutFreq ./Trajectory$(ChargeTemperSim ?  "_\${id}" : "").h5 position image create_group yes\n")
    end
    if (Restart)
        write(file, "dump_modify trajectory append yes\n")
    end

    
    write(file,"
    fix 1 all nve 
    fix 2 all langevin $(Temperature) $(Temperature) $(SimulationType=="Calvados2" || SimulationType=="ArashModell" ? 100 : 1000) $(rand(1:100_000)) ### Tstart, Tstop, Dampening in fs  (1ps dampening in paper), seed\n")
    if (Restart)
    write(file,"
    fix 3 all print $WriteOutFreq  \"\$(step) \$(time) \$(temp) \$(press) \$(etotal) \$(pe) \$(ke) \$(evdwl) \$(ecoul) \$(epair) \$(ebond) \$(eangle) \$(edihed)\" append \${energyfile} screen no\n")
    else
        write(file,"
    fix 3 all print $WriteOutFreq  \"\$(step) \$(time) \$(temp) \$(press) \$(etotal) \$(pe) \$(ke) \$(evdwl) \$(ecoul) \$(epair) \$(ebond) \$(eangle) \$(edihed)\" file \${energyfile} screen no\n")
    end

    write(file,"   
    restart 100000000 ./Restart/BackUp$(ChargeTemperSim ?  "_\${id}" : "").*
    timestep 10.0
    timer timeout 23:30:00")

    if ChargeTemperSim
        write(file,"   
    chargetemper $NTimeSteps $ChargeTemperSwapSteps \$q 2 $(rand(1:100_000)) $(rand(1:100_000)) ### timesteps swapfreq chargeinit temp_fix_id seed1 seed2 ")
    else
    write(file,"   
    run $(NTimeSteps) upto")
    end
    write(file,"
    write_restart ./sim$(ChargeTemperSim ?  "_\${id}" : "").restart")
    close(file)
end

function writeLammpsScript(fileName, AtomTypes; WriteOutFreq=100000)
    file = open(fileName, "w")

    write(file," ### Initialisation

    units real ### https://docs.lammps.org/stable/units.html
    dimension 3
    newton on ### might be faster if off depends on parallelization https://docs.lammps.org/stable/newton.html
    processors * * * grid numa #map zyx
    ### processors pre dimension, can optimize a lot here.
    boundary p p p # periodic in each dimension 
    atom_style full 
    bond_style harmonic/omp

    atom_modify id yes  #ids will be assigned to each atom 
    #atom_modify map hash 
    #atom_modify 10000 2 #resort ids every 10000 steps for cache improvements 


    # Reading file i easier
    read_data ./Start_conf.txt
    neighbor 2.0 multi
    neigh_modify every 5
    neigh_modify collection/interval 5 20 22 24 26 36


    ### Simulation Settings
    #Bonds
    bond_coeff 1 10.0 3.82   ### ????, k=10 kcal/mol/AA^2, sigma_0=3.82AA


    ### Coulomb and Ashbaugh Hatch Potential
    pair_style ah/cut/coul/cut 35 0.1 0.2 ### Dielectric constant
    dielectric 80.0 		### dielectric constant \n")

    for (id, ResType) in enumerate(AtomTypes)
        write(file, "pair_coeff $id $id $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType])  $(BioData.OneToHPSCharge[ResType])\n")
    end

    #=for (id, ResType) in enumerate(AtomTypes)
        for (id2, ResType2) in enumerate(AtomTypes)
            if id2<id
                continue
            end
            write(file, "pair_coeff $id $id2 ah $((BioData.OneToHPSUrrySigma[ResType]+BioData.OneToHPSUrrySigma[ResType2])/2.) $((BioData.OneToHPSUrryLambda[ResType]+BioData.OneToHPSUrryLambda[ResType])/2.) \n")
        end
        #write(file, "pair_coeff $id $id ah $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType]) \n")
    end=#
    write(file, "\npair_modify mix arithmetic shift yes\n")


    write(file, "run_style verlet

    angle_style none
    dihedral_style none
    dump trajectory2 all xyz $(WriteOutFreq) ./Trajectory.xyz 


    minimize 1.0e-4 1.0e-6 500 1000
    fix 1 all nve # temp 300.0 300.0 1000  ### 
    fix 2 all langevin 285 285 1000 $(rand(1:100_000))) ### Tstart, Tstop, Dampening in fs  (1ps dampening in paper), seed
    fix 3 all print $(WriteOutFreq)  \"\$(step) \$(time) \$(temp) \$(press) \$(etotal) \$(pe) \$(ke) \$(evdwl) \$(ecoul) \$(epair) \$(ebond) \$(eangle) \$(edihed)\" file ./energy.dat screen no

    restart $(WriteOutFreq) ./sim.restart
    thermo $(WriteOutFreq)
    timestep 10.0
    run 1000000000")
    close(file)
end
#Set start coordinates for the AA, by adding offset for x,y-position and offset in z-position if Seqindex%ProteinsPerLayer==0
function createStartingPosition(Sequences, BoxSize, LatticeSpacing=8.0, BondSpacing=3.8)
    NSeq = length(Sequences)
    MaxSeqLength =0 
    for Seq in Sequences
        MaxSeqLength = max(MaxSeqLength, length(Seq))
    end

    pos = zeros(NSeq, MaxSeqLength, 3)
    ProteinsPerLayer=ceil(sqrt(NSeq))
    NLayer = ceil(NSeq/ProteinsPerLayer)
    zoffset =(BoxSize[5]+BoxSize[6])/2. -((NLayer-1)/2.)*LatticeSpacing
    for (Seqindex, Seq) in enumerate(Sequences)
        xoffset = (BoxSize[1]+BoxSize[2])/2.
        xoffset += ( (Seqindex%ProteinsPerLayer) -(ProteinsPerLayer-1)/2.)*LatticeSpacing
        yoffset = (BoxSize[3]+BoxSize[4])/2.
        yoffset -= MaxSeqLength/2. *BondSpacing
        for (Aaid,_) in enumerate(Seq)
            pos[Seqindex, Aaid, 1] = xoffset
            pos[Seqindex, Aaid, 2] = yoffset
            pos[Seqindex, Aaid, 3] = zoffset
            yoffset += BondSpacing
        end
        if Seqindex%ProteinsPerLayer==0
            zoffset += LatticeSpacing
        end
    end
    return pos
end

function createDenseStartingPosition(Sequences, BoxSize, LatticeSpacing=8.0, BondSpacing=3.8)
    NSeq = length(Sequences)
    MaxSeqLength =0 
    for Seq in Sequences
        MaxSeqLength = max(MaxSeqLength, length(Seq))
    end

    pos = zeros(NSeq, MaxSeqLength, 3)
    xoffset= BoxSize[1]+0.5
    yoffset=BoxSize[3]+0.5
    zoffset =BoxSize[5]+0.5
    forward = 1.
    xforward = 1.
    ytmp= 0.
    for (Seqindex, Seq) in enumerate(Sequences)
        for (Aaid,_) in enumerate(Seq)
            ytmp =  yoffset+BondSpacing*forward
            if ytmp <(BoxSize[4] -0.5) && ytmp >(BoxSize[3] +0.5)
                yoffset = ytmp
            else
                xtmp =  xoffset+LatticeSpacing*xforward
                if xtmp<(BoxSize[2] -0.5) && xtmp >(BoxSize[1] +0.5)
                    xoffset = xtmp ### otherwise overlappting chains, non ideal harmonic bond will equilibrate fast
                else
                    xforward *=-1.
                    zoffset += LatticeSpacing
                    #xoffset = BoxSize[1]+0.5
                end
                forward *= -1.
            end
            pos[Seqindex, Aaid, 1] = xoffset
            pos[Seqindex, Aaid, 2] = yoffset
            pos[Seqindex, Aaid, 3] = zoffset
        end
    end
    return pos
end

function correctPositionInBounds(pos, boxSize, Sequences)
    box = boxSize/2.0
    box_neg = -1.0.*box

    for (SeqID, Seq) in enumerate(Sequences)
        for (atom,_) in enumerate(Seq)
            for i in 1:3
                if pos[SeqID, atom, i]>box[i]
                    pos[SeqID, atom, i] -= boxSize[i]
                end
                if pos[SeqID, atom, i]< box_neg[i]
                    pos[SeqID, atom, i] += boxSize[i]
                end
            end
        end
    end
    return pos
end

function rescalePositions(pos, boxSize, Sequences, image, rescale_factor=5.5/3.8)
    uw_pos = deepcopy(pos) 
    for i in 1:3
        uw_pos[:,:,i] .+= boxSize[i] .*image[:,:,i]
    end
    uw_pos ./= rescale_factor
    return correctPositionInBounds(uw_pos, boxSize, Sequences)
end

@doc raw"""
    DetermineYukawaInteractions(;SimulationType="", Temperature=300, SaltConcentration=-1)

Calculates the constants for the Yukawa potencial for given temperature and salt concentration.
    
**Arguments**
- `SimulationType::String`: Type of Simulation (e.g.: "Calvados2").

**Optional Arguments**:
- `Temperature::Float`: Temperature of the Simulation.
- `SaltConcentration::Float`: Salt concentration of the Simulation.

**Return**:
* `ϵ_r::Float`: Temperature-dependent dielectric constant.
* `κ::Float`: Invers of the Debye length.
"""
function DetermineYukawaInteractions(;SimulationType="", Temperature=300, SaltConcentration=-1)
    ### constants
	e = 1.602*10.0^-(19) ### Charge of electron
	e_0 = 8.854188*10.0^-(12) ### vacuum permitivity
	kb = 1.380*10.0^-(23) ### boltzmann constant
    NA = 6.02214086 * 10.0^23 # 1/mol Avogadro constant

    ### Saltconcentration in mol/L
    if SaltConcentration==-1 || !in(SimulationType, ["Calvados2","Calvados3"])
        κ = 10.0 # Å, Screening length
        ϵ_r = 80.0 # what ever? I think unitless, relative permitivity
    else
        ### taken from RESEARCH ARTICLE Improved predictions of phase behaviour of intrinsically disordered proteins by tuning the interaction range; Giulio Tesei , Kresten Lindorff-Larsen; 2022
        RT =  8.3145*Temperature*1e-3
	    ϵ_r= 5321.0/Temperature+233.76-0.9297*Temperature+1.417*10.0^-3*Temperature^2-8.292*10.0^-7*Temperature^3 ### relative permitivity
        #λ = e^2 /(4.0*pi*e_0*ϵ_r*kb*Temperature) ### m, Bjerrum length
        λ = 1.6021766^2/(4*pi*8.854188*ϵ_r)*6.022*1000/RT ### ?, Bjerrum length

        κ = sqrt(8.0*pi*λ*SaltConcentration*6.022/10)
        #κ = sqrt(1.0/(8.0*pi*λ*NA*SaltConcentration*1000.0))*10.0^10 ### Å Screening length, named D in the paper
    end
    return (ϵ_r, κ)
end


@doc raw"""
    CalvadosSetup(Sequences,AtomTypes,pH)

Define escential Parameters for the Simulation based on the Simulation Type.
    
**Arguments**
- `Sequences::Array{String}`: List of sequences of Proteins.
- `SimulationType::String`: Type of Simulation (e.g.: "Calvados2").
- `pH::Float`: pH-value of the system.

**Return**:
A tuple containing:
- `AtomTypes::Set{Char}`: Set of unique amino acid types in the provided sequences.
- `LongAtomTypes::Set{Char}`: Set of amino acid types where the first and last amino acid in a sequence are treated as distinct types.
- `AaToId::Dict{Char, Int32}`: Dictionary mapping each amino acid type to a unique ID number.
- `IdToAa::Dict{Int32, Char}`: Dictionary mapping each ID number to its corresponding amino acid type.
- `ResToLongAtomType::Dict{Tuple{Char, Bool}, Char}`: Dictionary mapping standard amino acids to their modified forms when at the beginning or end of a sequence.
- `LongAtomTypesToRes::Dict{Char, Tuple{Char, Bool}}`: Reverse mapping of `ResToLongAtomType`.
- `OneToCharge::Dict{Char, Float}`: Dictionary containing the charge values of amino acids, modified based on simulation type.
- `OneToMass::Dict{Char, Float}`: Dictionary containing the mass values of amino acids.
- `OneToSigma::Dict{Char, Float}`: Dictionary containing the sigma values of amino acids.
- `OneToLambda::Dict{Char, Float}`: Dictionary containing the lambda values of amino acids.
- `OneToHPSDihedral0110::Dict{Char, Any}`: Dictionary containing dihedral parameters for a specific configuration.
- `OneToHPSDihedral1001::Dict{Char, Any}`: Dictionary containing dihedral parameters for another configuration.

**Notes**:
- If `SimulationType` is "Calvados2", the first and last amino acids in each sequence are assigned different types to account for altered mass due to peptide bonding. Also the charge of histidine ('H') is adjusted based on the provided pH value using the formula: `1/(1+10^(pH-6))`.
"""
function CalvadosSetup(Sequences, AtomTypes, pH, AaToId, SimulationType)
    index_cnt = length(AtomTypes)
    LongAtomTypes=deepcopy(AtomTypes)
    ResToLongAtomType = Dict()
    NewSequences= deepcopy(Sequences)
    cnt = 96 ### use lower case ascii letters as modified residue types
    for (id,Sequence) in enumerate(Sequences)
        if( ~((Sequence[1], true) in  keys(ResToLongAtomType)))
            ResToLongAtomType[(Sequence[1], true)] = Char(cnt+=1) 
            push!(LongAtomTypes, ResToLongAtomType[(Sequence[1], true)] )#
            AaToId[ResToLongAtomType[(Sequence[1], true)] ] = index_cnt+=1
            NewSequences[id] =  ResToLongAtomType[(Sequence[1], true)]* Sequence[2:end]
        end
        if( ~((Sequence[end], false) in  keys(ResToLongAtomType)))
            ResToLongAtomType[(Sequence[end], false)] = Char(cnt+=1)
            push!(LongAtomTypes, ResToLongAtomType[(Sequence[end], false)])
            AaToId[ResToLongAtomType[(Sequence[end], false)] ] = index_cnt+=1
            NewSequences[id] = NewSequences[id][1:end-1]* ResToLongAtomType[(Sequence[end], false)]
        end
    end

    Sequences.=NewSequences
    LongAtomTypes = union(LongAtomTypes, AtomTypes)
    LongAtomTypesToRes=Dict( (v => k) for (k, v) in ResToLongAtomType)

    OneToCharge = deepcopy(BioData.OneToHPSCharge)
    OneToMass = deepcopy(BioData.AaToWeight)
    OneToSigma  = deepcopy(BioData.OneToHPSCalvadosSigma)
    OneToHPSDihedral0110 = deepcopy(BioData.OneToHPSDihedral0110)
    OneToHPSDihedral1001 = deepcopy(BioData.OneToHPSDihedral1001)
    OneToCharge['H'] = 1.0/(1+10^(pH-6))  ### HIS is pH-Dependent for Calvados2/Calvados3

    if SimulationType=="Calvados2"
        OneToLambda = deepcopy(BioData.OneToCalvados2Lambda)
    elseif SimulationType=="Calvados3"
        OneToLambda = deepcopy(BioData.OneToCalvados3Lambda)
    end

    for e in LongAtomTypes ### added the charge modifications for the first/last amino acids
        if ~(e in keys(OneToCharge))
            (AA, front) = LongAtomTypesToRes[e]
            OneToCharge[e] = front ? OneToCharge[AA] +1 :  OneToCharge[AA] -1
            OneToMass[e] = front ? OneToMass[AA] +2.0 :  OneToMass[AA] +16
            OneToSigma[e] = OneToSigma[AA]
            OneToLambda[e] = OneToLambda[AA]
            OneToHPSDihedral0110[e] = OneToHPSDihedral0110[AA]
            OneToHPSDihedral1001[e] = OneToHPSDihedral1001[AA]
        end
    end
    return (AtomTypes, LongAtomTypes, AaToId,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)
end

#=
@doc raw"""
    DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma)

Define escential Parameters for the Simulation based on the Simulation Type.
    
**Arguments**
- `Sequences::Array{String}`: List of sequences of Proteins.
- `SimulationType::String`: Type of Simulation (e.g.: "Calvados2").
- `pH::Float`: pH-value of the system.

**Return**:
A tuple containing:
- `AtomTypes::Set{Char}`: Set of unique amino acid types in the provided sequences.
- `LongAtomTypes::Set{Char}`: Set of amino acid types where the first and last amino acid in a sequence are treated as distinct types.
- `AaToId::Dict{Char, Int32}`: Dictionary mapping each amino acid type to a unique ID number.
- `IdToAa::Dict{Int32, Char}`: Dictionary mapping each ID number to its corresponding amino acid type.
- `ResToLongAtomType::Dict{Tuple{Char, Bool}, Char}`: Dictionary mapping standard amino acids to their modified forms when at the beginning or end of a sequence.
- `LongAtomTypesToRes::Dict{Char, Tuple{Char, Bool}}`: Reverse mapping of `ResToLongAtomType`.
- `OneToCharge::Dict{Char, Float}`: Dictionary containing the charge values of amino acids, modified based on simulation type.
- `OneToMass::Dict{Char, Float}`: Dictionary containing the mass values of amino acids.
- `OneToSigma::Dict{Char, Float}`: Dictionary containing the sigma values of amino acids.
- `OneToLambda::Dict{Char, Float}`: Dictionary containing the lambda values of amino acids.
- `OneToHPSDihedral0110::Dict{Char, Any}`: Dictionary containing dihedral parameters for a specific configuration.
- `OneToHPSDihedral1001::Dict{Char, Any}`: Dictionary containing dihedral parameters for another configuration.

**Notes**:
- If `SimulationType` is "Calvados2", the first and last amino acids in each sequence are assigned different types to account for altered mass due to peptide bonding. Also the charge of histidine ('H') is adjusted based on the provided pH value using the formula: `1/(1+10^(pH-6))`.
"""
function CalvadosSetup(Sequences, AtomTypes, pH, AaToId, SimulationType)
    index_cnt = length(AtomTypes)
    LongAtomTypes=deepcopy(AtomTypes)
    ResToLongAtomType = Dict()
    NewSequences= deepcopy(Sequences)
    cnt = 96 ### use lower case ascii letters as modified residue types
    for (id,Sequence) in enumerate(Sequences)
        if( ~((Sequence[1], true) in  keys(ResToLongAtomType)))
            ResToLongAtomType[(Sequence[1], true)] = Char(cnt+=1) 
            push!(LongAtomTypes, ResToLongAtomType[(Sequence[1], true)] )#
            AaToId[ResToLongAtomType[(Sequence[1], true)] ] = index_cnt+=1
            NewSequences[id] =  ResToLongAtomType[(Sequence[1], true)]* Sequence[2:end]
        end
        if( ~((Sequence[end], false) in  keys(ResToLongAtomType)))
            ResToLongAtomType[(Sequence[end], false)] = Char(cnt+=1)
            push!(LongAtomTypes, ResToLongAtomType[(Sequence[end], false)])
            AaToId[ResToLongAtomType[(Sequence[end], false)] ] = index_cnt+=1
            NewSequences[id] = NewSequences[id][1:end-1]* ResToLongAtomType[(Sequence[end], false)]
        end
    end

    Sequences.=NewSequences
    LongAtomTypes = union(LongAtomTypes, AtomTypes)
    LongAtomTypesToRes=Dict( (v => k) for (k, v) in ResToLongAtomType)

    OneToCharge = deepcopy(BioData.OneToHPSCharge)
    OneToMass = deepcopy(BioData.AaToWeight)
    OneToSigma  = deepcopy(BioData.OneToHPSCalvadosSigma)
    OneToHPSDihedral0110 = deepcopy(BioData.OneToHPSDihedral0110)
    OneToHPSDihedral1001 = deepcopy(BioData.OneToHPSDihedral1001)
    OneToCharge['H'] = 1.0/(1+10^(pH-6))  ### HIS is pH-Dependent for Calvados2/Calvados3

    if SimulationType=="Calvados2"
        OneToLambda = deepcopy(BioData.OneToCalvados2Lambda)
    elseif SimulationType=="Calvados3"
        OneToLambda = deepcopy(BioData.OneToCalvados3Lambda)
    end

    for e in LongAtomTypes ### added the charge modifications for the first/last amino acids
        if ~(e in keys(OneToCharge))
            (AA, front) = LongAtomTypesToRes[e]
            OneToCharge[e] = front ? OneToCharge[AA] +1 :  OneToCharge[AA] -1
            OneToMass[e] = front ? OneToMass[AA] +2.0 :  OneToMass[AA] +16
            OneToSigma[e] = OneToSigma[AA]
            OneToLambda[e] = OneToLambda[AA]
            OneToHPSDihedral0110[e] = OneToHPSDihedral0110[AA]
            OneToHPSDihedral1001[e] = OneToHPSDihedral1001[AA]
        end
    end
    return (AtomTypes, LongAtomTypes, AaToId,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)
end
=#

@doc raw"""
    DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma)

Defines essential Parameters for the simulation based on the type of the simulation.
    
**Arguments**
- `Sequences::Array{String}`: List of sequences of proteins.
- `SimulationType::String`: Type of simulation (e.g.: "Calvados2").
- `pH::Float`: pH-value of the system.

**Optional Arguments**:
- `OneToChargeDef::Dict()`: Dictionary defining the charge for the Aminoacids.
- `OneToLambdaDef::Dict()`: Dictionary defining the Lambda for the Aminoacids.
- `OneToSigmaDef::Dict()`: Dictionary defining the Sigma for the Aminoacids.

**Return**:
A tuple containing:
- `AtomTypes::Set{Char}`: Set of unique amino acid types in the provided sequences.
- `LongAtomTypes::Set{Char}`: Set of amino acid types where the first and last amino acid in a sequence are treated as distinct types.
- `AaToId::Dict{Char, Int32}`: Dictionary mapping each amino acid type to a unique ID number.
- `IdToAa::Dict{Int32, Char}`: Dictionary mapping each ID number to its corresponding amino acid type.
- `ResToLongAtomType::Dict{Tuple{Char, Bool}, Char}`: Dictionary mapping standard amino acids to their modified forms when at the beginning or end of a sequence.
- `LongAtomTypesToRes::Dict{Char, Tuple{Char, Bool}}`: Reverse mapping of `ResToLongAtomType`.
- `OneToCharge::Dict{Char, Float}`: Dictionary containing the charge values of amino acids, modified based on simulation type.
- `OneToMass::Dict{Char, Float}`: Dictionary containing the mass values of amino acids.
- `OneToSigma::Dict{Char, Float}`: Dictionary containing the sigma values of amino acids.
- `OneToLambda::Dict{Char, Float}`: Dictionary containing the lambda values of amino acids.
- `OneToHPSDihedral0110::Dict{Char, Any}`: Dictionary containing dihedral parameters for a specific configuration.
- `OneToHPSDihedral1001::Dict{Char, Any}`: Dictionary containing dihedral parameters for another configuration.

**Notes**:
- If `SimulationType` is "Calvados2", the first and last amino acids in each sequence are assigned different types to account for altered mass due to peptide bonding. Also the charge of histidine ('H') is adjusted based on the provided pH value using the formula: `1/(1+10^(pH-6))`.
- If `SimulationType` is "HPS-Alpha", predefined values for charge, lambda, and sigma are used.
- If an unknown `SimulationType` is provided, the function falls back on the user-supplied dictionaries for charge, lambda, and sigma values.
"""
function DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma)
    #Define Aminoacids to ID Dict
    AtomTypes= Set(join(Sequences))
    AaToId = Dict{Char,Int32}()
    for (index, value) in enumerate(AtomTypes)
        AaToId[value]=index
    end

    LongAtomTypes=deepcopy(AtomTypes)
    ResToLongAtomType = Dict()
    LongAtomTypesToRes=Dict()
    OneToMass = deepcopy(BioData.AaToWeight)
    OneToHPSDihedral0110 = deepcopy(BioData.OneToHPSDihedral0110)
    OneToHPSDihedral1001 = deepcopy(BioData.OneToHPSDihedral1001)
    if SimulationType=="HPS-Alpha"
        OneToCharge = deepcopy(BioData.OneToHPSCharge)
        OneToLambda = deepcopy(BioData.OneToHPSUrryLambda)
        OneToSigma  = deepcopy(BioData.OneToHPSCalvadosSigma)
    elseif SimulationType=="Calvados2"||SimulationType=="Calvados3"
        (AtomTypes, LongAtomTypes, AaToId,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)=CalvadosSetup(Sequences, AtomTypes, pH, AaToId, SimulationType)
    else ### take the ones which are supplied
        OneToCharge = deepcopy(OneToChargeDef)
        OneToLambda = deepcopy(OneToLambdaDef)
        OneToSigma = deepcopy(OneToSigmaDef)
    end
    IdToAa=Dict( (v => k) for (k, v) in AaToId)

    return (AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)
end

function getImageCopyNumber(pos, boxSize, Sequences)
    diff = zeros(eltype(pos), 3)
    image = zeros(Int32, size(pos))
    tmp_pos = deepcopy(pos)

    for (SeqID, Seq) in enumerate(Sequences)
        for (atom,_) in enumerate(Seq)
            if atom ==1 continue end
            for i in 1:3
                diff[i] = tmp_pos[SeqID, atom, i]-tmp_pos[SeqID, atom-1, i]

                if abs(diff[i]-boxSize[i])<abs(diff[i]) # || abs(diff[i]+boxSize[i])<abs(diff[i])
                    image[SeqID, atom,i] =  image[SeqID, atom-1,i]-1 #-1*convert(Int32,(round((diff[i]/(boxSize[i]/2.)/2.0))))  + image[SeqID, atom-1,i]
                    #tmp_pos[SeqID, atom, i] +=
                elseif abs(diff[i]+boxSize[i])<abs(diff[i])
                    image[SeqID, atom,i] =  image[SeqID, atom-1,i]+1
                    #image[SeqID, atom,i] = -1*convert(Int32,(round((diff[i]/(boxSize[i]/2.)/2.0))))  + image[SeqID, atom-1,i]
                else
                    image[SeqID, atom,i] = image[SeqID, atom-1,i]
                end
            end
        end
    end

    return image
end

function startConfigurationSetup(Sequences,SimulationType,pH,OneToChargeDef,OneToLambdaDef,OneToSigmaDef,MixingRule="1-1001-1")
    #Define Number of all Aminoacids, Bonds, Angles and Dihedrals, set AlphaAddition
    NAtoms=0
    NBonds=0
    NAngles=0
    NDihedrals=0
    for Seq in Sequences
        NAtoms += length(Seq)
        NBonds += length(Seq)-1
        NAngles += length(Seq)-2
        NDihedrals += length(Seq)-3
    end

    AlphaAddition=false
    if SimulationType=="Calvados2+Alpha"
        AlphaAddition = true
        SimulationType="Calvados2"
    elseif SimulationType=="Calvados3+Alpha"
        AlphaAddition = true
        SimulationType="Calvados3"
    else
        AlphaAddition=false
    end

    (AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001) =  DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=OneToChargeDef, OneToLambdaDef=OneToLambdaDef, OneToSigmaDef=OneToSigmaDef)
    #Define length of all chains
    NAtomTypes = length(LongAtomTypes)

    #if AlphaAddition then determine the Dihedral
    dihedral_short_map=Dict()
    dihedral_list = zeros(Int32, (0,0))
    dihedral_long_map=Dict()
    dihedral_eps=[]
    if AlphaAddition
        (dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list) = determineDihedrals(Sequences, AtomTypes, AaToId, OneToHPSDihedral0110, OneToHPSDihedral1001, MixingRule)
        NDihedralsTypes = length(dihedral_eps)
    end
    return NAtoms,NBonds,NAngles,NDihedrals,AlphaAddition,SimulationType,AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, NAtomTypes, dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list
end
@doc raw"""
    writeHOOMD(Sequences,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,Device,yk_cut,ah_cut,pH,domain,NAtoms,NBonds,NAngles,NDihedrals,dihedral_short_map,dihedral_list, ENM)

Writes the datas for the simulation in different files with HOOMD.
    
**Arguments**
- from writeStartConfiguration()

**Creates**:
* Writes the datas for the simulation.
"""
function writeHOOMD(BasePath, Sequences,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,Device,yk_cut,ah_cut,pH,domain,NAtoms,NBonds,NAngles,NDihedrals,dihedral_short_map,dihedral_list, ENM)
        mkpath("$(BasePath)/HOOMD_Setup")
        WriteHOOMDSequences("$(BasePath)/HOOMD_Setup/Sequences.txt", Sequences)
        WriteHOOMDParticlesInput("$(BasePath)/HOOMD_Setup/Particles.txt", pos,  OneToCharge, AaToId,Sequences, OneToMass, OneToSigma, image)
        if AlphaAddition
            WriteDihedrals("$(BasePath)/HOOMD_Setup/DihedralMap.txt", dihedral_long_map, dihedral_eps)
        end
        (ϵ_r, κ) = DetermineYukawaInteractions(;SimulationType=SimulationType, Temperature=Temperature, SaltConcentration=SaltConcentration)

        BoxLength = [BoxSize[2]-BoxSize[1],BoxSize[4]-BoxSize[3],BoxSize[6]-BoxSize[5] ]
        WriteParams("$(BasePath)/HOOMD_Setup/Params.txt", StartFileName, Temperature, NSteps, WriteOutFreq, 0.01, BoxLength/10.0, rand(1:65535), UseAngles=AlphaAddition;ϵ_r=ϵ_r, κ=κ,Device=Device, yk_cut=yk_cut/10.0, ah_cut=ah_cut/10.0, ionic=SaltConcentration, pH=pH, SimType=SimulationType, domain=domain,Create_Start_Config=true) ### BoxLength has to be convert to nm
        WriteDictionaries("$(BasePath)/HOOMD_Setup/Dictionaries.txt", OneToCharge, AaToId, OneToMass, OneToSigma, OneToLambda)
        WriteENM_HOOMD_Indices("$(BasePath)/HOOMD_Setup/ENM_indices.txt", ENM)
        InputMasses = [OneToMass[res] for res in join(Sequences)]
        InputCharges = [OneToCharge[res] for res in join(Sequences)]
        writeGSDStartFile("$BasePath$StartFileName.gsd", NAtoms, NBonds, NAngles, NDihedrals,BoxLength, pos, AaToId,Sequences,image, InputMasses, InputCharges, dihedral_short_map, dihedral_list, OneToSigma, AlphaAddition, SimulationType, domain, ENM)    
end

@doc raw"""
    writeHPSLammps(fileName,Sequences,AtomTypes,LongAtomTypes,LongAtomTypesToRes,InitStyle,ChargeTemperSteps,ChargeTemperSwapSteps,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,pH,NAtoms,NBonds,NAngles,NDihedrals,Info,NAtomTypes)

Writes the datas for the simulation in different files.
    
**Arguments**
- from writeStartConfiguration()

**Creates**:
* Writes the datas for the simulation.
"""
function writeHPSLammps(fileName,Sequences,AtomTypes,LongAtomTypes,LongAtomTypesToRes,InitStyle,ChargeTemperSteps,ChargeTemperSwapSteps,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,pH,NAtoms,NBonds,NAngles,NDihedrals,Info,NAtomTypes)
    writeHPSLammpsScript( fileName*".lmp",StartFileName, AtomTypes, LongAtomTypes, AaToId, LongAtomTypesToRes, OneToCharge, OneToSigma, OneToLambda, dihedral_eps, InitStyle, SimulationType, Temperature, AlphaAddition, false, NSteps; SaltConcentration=SaltConcentration, pH=pH, ChargeTemperSteps=ChargeTemperSteps, ChargeTemperSwapSteps=ChargeTemperSwapSteps,WriteOutFreq=WriteOutFreq)
    writeHPSLammpsScript( fileName*"_restart.lmp",StartFileName, AtomTypes, LongAtomTypes, AaToId, LongAtomTypesToRes, OneToCharge, OneToSigma, OneToLambda, dihedral_eps, InitStyle, SimulationType, Temperature, AlphaAddition, true, NSteps; SaltConcentration=SaltConcentration, pH=pH, ChargeTemperSteps=ChargeTemperSteps, ChargeTemperSwapSteps=ChargeTemperSwapSteps,WriteOutFreq=WriteOutFreq)

    file = open(StartFileName*".txt", "w");
    write(file, Info)

    write(file, "\n\t $NAtoms \t atoms\n")
    write(file, "\t $NBonds \t bonds\n")
    write(file, "\t $NAngles \t angles\n")
    write(file, "\t $NDihedrals \t dihedrals\n\n")

    write(file, "\t $NAtomTypes \t atom types\n")
    write(file, "\t 1 \t bond types\n")
    write(file, "\t 1 \t angle types\n")
    write(file, "\t $NDihedralsTypes \t dihedral types\n\n")

    write(file,"\t $(BoxSize[1]) \t $(BoxSize[2]) \t xlo xhi\n")
    write(file,"\t $(BoxSize[3]) \t $(BoxSize[4]) \t ylo yhi\n")
    write(file,"\t $(BoxSize[5]) \t $(BoxSize[6]) \t zlo zhi\n\n")

    write(file, "Masses \n#\n")
    for (index, value) in enumerate(LongAtomTypes)
        write(file, "\t $(AaToId[value]) \t $(OneToMass[value]) \t ### $(value in AtomTypes ? value : LongAtomTypesToRes[value] ) \n")
    end

    write(file,"\nAtoms\n# A comment that is needed to read stuff accurately\n")

    atomid=0
    moleculeid=0
    for (SeqId, Sequence) in enumerate(Sequences)
        moleculeid+=1
        for (ResId,Res) in enumerate(Sequence)
            atomid +=1
            write(file, "\t $atomid \t $(@sprintf("%6i", moleculeid)) \t  $(@sprintf("%2i", AaToId[Res])) \t  $(@sprintf("%3.5f", OneToCharge[Res])) \t $(@sprintf("%5.5f", pos[SeqId,ResId,1])) \t $(@sprintf("%5.5f", pos[SeqId,ResId,2])) \t $(@sprintf("%5.5f", pos[SeqId,ResId,3])) \t $(@sprintf("%i", image[SeqId,ResId,1])) \t $(@sprintf("%i", image[SeqId,ResId,2])) \t $(@sprintf("%i", image[SeqId,ResId,3]))\n")
        end
    end 

    write(file,"\nBonds\n#\n")
    bonds = getBonds(Sequences;M=2)
    for bid in axes(bonds,1)
        write(file, "\t $bid \t 1 \t  $(bonds[bid, 1]) \t  $(bonds[bid, 1]) \n")
    end
    write(file, "\n\n")


    write(file,"\nAngles\n#\n")

    angles = getBonds(Sequences;M=3)
    for bid in axes(bonds,1)
        write(file, "\t $bid \t 1 \t  $(angles[bid, 1]) \t  $(angles[bid, 1]) \t  $(angles[bid, 2])\n")
    end
    write(file, "\n\n")


    write(file,"\nDihedrals\n#\n")
    atomid=0
    dihedralid = 0
    for (SeqId, Seq) in enumerate(Sequences)
        for (ResId,Res) in enumerate(Seq)
            atomid +=1
            if (ResId>(length(Seq)-3) )
                continue
            else
                dihedralid += 1
                if MixingRule=="1-1001-1"
                    Res_min = (ResId-1)>=1 ? AaToId[Seq[ResId-1]] : 0
                    Res_max = (ResId)<=(length(Seq)-4) ? AaToId[Seq[ResId+4]] : -1
                    key = (sort([Res_min,AaToId[Res], AaToId[Seq[ResId+3]], Res_max]))
                elseif MixingRule=="1001"
                    key = (sort([AaToId[Res], AaToId[Seq[ResId+3]]]))
                elseif MixingRule=="0110"
                    key = (sort([AaToId[Seq[ResId+1]], AaToId[Seq[ResId+2]]]))
                end
                write(file, "\t $dihedralid \t $(dihedral_long_map[key]) \t  $atomid \t  $(atomid+1) \t $(atomid+2) \t $(atomid+3)\n")
            end
        end
    end
    close(file)    
end

@doc raw"""
    writeStartConfiguration(fileName, StartFileName, Info, Sequences, BoxSize,NSteps=100_000_000; SimulationType="Calvados2", Temperature=300,MixingRule="1-1001-1", Pos =zeros(Float32, 0),InitStyle="Slab", SaltConcentration=0.15, pH=6.0, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, HOOMD=false, OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma,WriteOutFreq=100_000, Device="GPU", yk_cut=40.0, ah_cut=20.0)

Writes the start configuration for a molecular dynamics simulation.
    
**Arguments**
- `fileName::String`: Name of the output file.
- `StartFileName::String`: Name of the initial configuration file.
- `Info::String`: Additional information about the simulation.
- `Sequences::Vector{String}`: List of amino acid sequences corresponding to the proteins.
- `BoxSize::Vector{Float}`: A vector of minmal/maximal box dimensions in each axis. ([x_min, x_max, y_min, y_max, z_min, z_max]).
- `NSteps::Int`: Number of simulation steps (default: 100,000,000).
- `SimulationType::String`: IDP model used for simulations.(default: "Calvados2").
- `Temperature::Float`: Temperature in Kelvin (default: 300).
- `MixingRule::String`: Mixing rule for optional dihedral potential for Calvados models.
- `Pos::Vector{Float}`: Initial positions. (default: empty array).
- `InitStyle::String`: Initialization style, e.g., "Slab" or "Pos".
- `SaltConcentration::Float`: Salt concentration in M (default: 0.15).
- `pH::Float`: pH level (default: 6).
- `ChargeTemperSteps`: List of charge tempering steps.
- `ChargeTemperSwapSteps::Int`: Swap steps for charge tempering.
- `HOOMD::Boolean`: Boolean to enable HOOMD compatibility.
- `OneToChargeDef::Dict`: Dictionary defining the amount of charge for each one letter atom type.
- `OneToSigmaDef::Dict`: Dictionary defining the sigma parameter in LJ potentials for each one letter atom type.
- `OneToLambdaDef::Dict`: Dictionary defining the lambda paramter in Ashbaugh-Hatch potentials for each one letter atom type.
- `WriteOutFreq::Int`: Frequency of writing output (default: 100,000).
- `Device::String`: Computational device, "CPU"/"GPU" (default: "GPU").
- `yk_cut::Float`: Cutoff distances for yukawa potential.
- `ah_cut::Float`: Cutoff distances for ashbaugh hatch potential.
- `domain::List(Int)`: Domains in which the ENM is active.
- `ENM::Tuple`: Data that are nessesary for ENM (bond name, id, group).

**Creates**:
* Writes data files with the start configuration.
"""
function writeStartConfiguration(BasePath, fileName, StartFileName, Info, Sequences, BoxSize,NSteps=100_000_000; SimulationType="Calvados2", Temperature=300,MixingRule="1-1001-1", Pos =zeros(Float32, 0),InitStyle="Slab", SaltConcentration=0.15, pH=6, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, HOOMD=false, OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma,WriteOutFreq=100_000, Device="GPU", yk_cut=40.0, ah_cut=20.0, domain=Array([]), ENM)

    ChargeTemperSim=length(ChargeTemperSteps)>0

    NAtoms,NBonds,NAngles,NDihedrals,AlphaAddition,SimulationType,AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, NAtomTypes, dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list = startConfigurationSetup(Sequences,SimulationType,pH,OneToChargeDef,OneToLambdaDef,OneToSigmaDef,MixingRule)


    #Set start coordinates for the AA, with different variations
    if InitStyle=="Slab"
        #By adding offset for x,y-position and offset in z-position if Seqindex%ProteinsPerLayer==0
        pos = createStartingPosition(Sequences, BoxSize)
    elseif (InitStyle=="Pos")
        #By getting Initial Positions
        cnt = 1
        pos_init = Pos

        pos = zeros(length(Sequences), maximum(length.(Sequences)), 3)
        for (Seqindex, Seq) in enumerate(Sequences)
            for (Aaid,_) in enumerate(Seq)
                pos[Seqindex, Aaid, 1] = pos_init[cnt,1]
                pos[Seqindex, Aaid, 2] = pos_init[cnt,2]
                pos[Seqindex, Aaid, 3] = pos_init[cnt,3]
                cnt+=1
            end
        end
    else
        #By getting Initial Positions
        pos = createDenseStartingPosition(Sequences, BoxSize)
    end

    #AltBox=[BoxLengthShort,BoxLengthLong,BoxLengthShort]
    AltBox = [BoxSize[2]-BoxSize[1], BoxSize[4]-BoxSize[3], BoxSize[6]-BoxSize[5]]
    #Periodic Boundery, if outside -> define coordinates inside and set image
    pos = Setup.correctPositionInBounds(pos, AltBox, Sequences) ### poly coordinates are not in box sometimes
    pos = Setup.correctPositionInBounds(pos, AltBox, Sequences) ### 2 times is the charm....
    pos = Setup.correctPositionInBounds(pos, AltBox, Sequences) ### 3 times is the charm....

    image = Setup.getImageCopyNumber(pos, AltBox, Sequences)

    domain = []
    #Write all Inputs, Parameters (Yukawa Interaction with Debye-Hückle), Dictionaries and the Start-File
    if HOOMD
        writeHOOMD(BasePath, Sequences,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,Device,yk_cut,ah_cut,pH,domain,NAtoms,NBonds,NAngles,NDihedrals,dihedral_short_map,dihedral_list, ENM)
    else
        writeHPSLammps(fileName,Sequences,AtomTypes,LongAtomTypes,LongAtomTypesToRes,InitStyle,ChargeTemperSteps,ChargeTemperSwapSteps,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,pH,NAtoms,NBonds,NAngles,NDihedrals,Info,NAtomTypes)
    end
end

function getBonds(Sequences::Vector{String}; M=2)
    N = sum(length.(Sequences).-(M-1))
    bonds = zeros(Int32, (N, M))
    bondid = 0
    atomid = 0
    for (SeqId, Sequence) in enumerate(Sequences)
        for (ResId,Res) in enumerate(Sequence)
            if ResId>(length(Sequence)-M+1)
                atomid += 1
                continue
            else
                bondid += 1
                for (id, i) in enumerate(0:M-1)
                    bonds[bondid,id] = atomid + i 
                end
                atomid += 1
            end
        end
    end 
    return bonds
end

function writeSlurmScript(fileName, Proteins, Path, RunsPerProtein, OMPCores, Restart=false)
    file = open(fileName, "w");
    write(file,"#!/bin/bash\nexport lmp=/p/home/jusers/witzky1/juwels/LAMMPS_HREMD/bin/lmp\n\n")
    for (ProtId, Protein) in enumerate(Proteins)
        for Run in 1:RunsPerProtein[ProtId]
            write(file,"cd ./$(Protein)/RUN_$( lpad(Run,3,"0"))/\n")
            if (Restart)
                write(file,"srun -n 1 --exclusive --exact --cpus-per-task=$(OMPCores[ProtId]) \$lmp -in ./$(Protein)_restart.lmp -sf omp -package omp $(OMPCores[ProtId]) &\n")
            else
                write(file,"srun -n 1 --exclusive --exact --cpus-per-task=$(OMPCores[ProtId]) \$lmp -in ./$(Protein).lmp -sf omp -package omp $(OMPCores[ProtId]) & \n")
            end
            write(file, "cd ../../\n\n")
        end
    end
    write(file, "wait ")
    close(file)
end

function writeCollectedSlurmScript(Path, Proteins, RelPaths,MPICores,OMPCores; ProteinsPerSlurmFile=2, Restarts=6, SlurmAccName="rsproteins", LmpAddOn="", CoresPerNodes=48, JobName="RS")
    Nodes=Int32(ceil(maximum(MPICores)*ProteinsPerSlurmFile/CoresPerNodes))
    slurm_file = open(Path*"run_all.sh", "w");
    write(slurm_file,"#!/bin/bash\nACCOUNT=$(SlurmAccName)\nRUNTIME=\"24:00:00\"\nPARTITION=\"batch\"\nNAME=\"$JobName\"\nMAIL=\" --mail-type=END,FAIL --mail-user=ywitzky@students.uni-mainz.de\"\nmodule load Intel\nmodule load IntelMPI\nmodule load HDF5/1.14.2-serial\nexport BASEPATH=\$PWD\n\n\n")
    cnt=0
    StartFile = open(Path*"test.sh", "w")
    RestartFile = open(Path*"test.sh", "w")
    for (RunNum, Protein) in enumerate(Proteins)
        if(mod(RunNum,ProteinsPerSlurmFile)==1  || ProteinsPerSlurmFile==1)
            if(RunNum÷ProteinsPerSlurmFile>0)
                write(StartFile, "\n\nwait ")
                write(RestartFile, "\n\nwait ")
                close(StartFile)
                close(RestartFile)
                if ProteinsPerSlurmFile>1
                    write(slurm_file,"\njid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL StartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
                    for ind in 1:Restarts
                        write(slurm_file,"jid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL --dependency=afterok:\${jid$(cnt-1)##* } RestartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
                    end
                end
            end
            StartFile = open(Path*"StartScript_$(RunNum÷ProteinsPerSlurmFile)", "w");
            write(StartFile,"#!/bin/bash\nexport lmp=/p/home/jusers/witzky1/juwels/LAMMPS_HREMD/bin/lmp\n\n")
            RestartFile = open(Path*"RestartScript_$(RunNum÷ProteinsPerSlurmFile)", "w");
            write(RestartFile,"#!/bin/bash\nexport lmp=/p/home/jusers/witzky1/juwels/LAMMPS_HREMD/bin/lmp\n\n")
        end
        write(StartFile,"cd ../$(RelPaths[RunNum]) \n")
        write(RestartFile,"cd ../$(RelPaths[RunNum]) \n")\
        write(RestartFile,"srun -n $(MPICores[RunNum]) --exclusive --exact --cpus-per-task=$(OMPCores[RunNum]) \$lmp -in \$BASEPATH/../$(RelPaths[RunNum])/$(Protein)_restart.lmp $(LmpAddOn)  -sf omp -package omp $(OMPCores[RunNum]) &\n")
        write(StartFile,"srun -n $(MPICores[RunNum]) --exclusive --exact --cpus-per-task=$(OMPCores[RunNum]) \$lmp -in \$BASEPATH/../$(RelPaths[RunNum])/$(Protein).lmp $(LmpAddOn)  -sf omp -package omp $(OMPCores[RunNum]) & \n")
        write(StartFile, "cd \$BASEPATH\n\n")
        write(RestartFile, "cd \$BASEPATH\n\n")
    end

    if(length(Proteins)÷ProteinsPerSlurmFile>0 || ProteinsPerSlurmFile==1 )

        RunNum=length(Proteins)+1
        write(StartFile, "\n\nwait ")
        write(RestartFile, "\n\nwait ")
        close(StartFile)
        close(RestartFile)
        write(slurm_file,"\njid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL StartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
        for ind in 1:Restarts
            write(slurm_file,"jid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL --dependency=afterok:\${jid$(cnt-1)##* } RestartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
        end
        #end
    end
    
    close(slurm_file)
end

@doc raw"""
    BuildENMModel(Sim::HPSAnalysis.SimData{T,I}, DomainDict, Proteins, Sequences, ProteinJSON) where {T<:Real, I<:Integer}

Calculate the Indices, that are nessesary to creat a start file im HOOMD.
    
**Arguments**
- `Sim::HPSAnalysis.SimData{T,I}`: The simulation datas.
- `DomainDict`: The Domains in which the ENM is active.
- `Proteins`: List of Protein Names.
- `Sequences`: The Sequences of the Proteins.
- `ProteinJSON`: AlphaFold data of the Proteins.

**Return**:
- Number of bonds
- Vector of bond type names 
- ID connecting bond tuple to bond type
- Vector of tuples defining all bonds
- Dict{String, Dict{Symbol, Float64}} which defines the bonds as used in HOOMD.
"""
function BuildENMModel(Sim::HPSAnalysis.SimData{T,I}, DomainDict, Proteins, Sequences, ProteinJSON) where {T<:Real, I<:Integer} 

    ConstraintDict, Backbone_correction_Dict = DetermineCalvados3ENMfromAlphaFold(Sim.BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0)
    HOOMD_Indices = ComputeHOOMD_ENMIndices(ConstraintDict, Backbone_correction_Dict, Sequences, Proteins)
    UnfoldedRegions =  GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)
    All_Indices = CombineBackboneAndENM(Proteins, Sequences, HOOMD_Indices, UnfoldedRegions, Backbone_correction_Dict)

    return All_Indices
end


@doc raw"""
    CombineBackboneAndENM(Proteins, Sequences, HOOMD_Indices, UnfoldedRegions, BackboneCorrectionDict)

Combines the bond length potential of backbone and ENM.
    
**Arguments**
- `Proteins`: List of Protein Names.
- `Sequences`: The Sequences of the Proteins.
- `HOOMD_Indices`: Tuple containing the info of ENM as defined by this functions return.
- `UnfoldedRegions`: Dictionary defining the unfolded regions.
- `BackboneCorrectionDict`::Dictionary with atoms, bond length and id which will be connected via backbone.

**Return**:
- Number of bonds
- Vector of bond type names 
- ID connecting bond tuple to bond type
- Vector of tuples defining all bonds
- Dict{String, Dict{Symbol, Float64}} which defines the bonds as used in HOOMD.
"""
function CombineBackboneAndENM(Proteins, Sequences, HOOMD_Indices, UnfoldedRegions, BackboneCorrectionDict)
    offsets = vcat([0],cumsum(length.(Sequences)))

    BackboneDiffLengths = Dict([prot=> [(i,j) for (i,j,_) in BackboneCorrectionDict[prot]] for prot in Set(Proteins)]) 

    N_Bonds = sum([sum([b-a for (a,b) in UnfoldedRegions[prot]]) for prot in Proteins])
    BB_ID = zeros(Int32, N_Bonds)
    BB_groups = Vector{Tuple{Int32, Int32}}()


    (ENM_Bonds, ENM_types, ENM_typeid, ENM_groups, harmonic) = HOOMD_Indices
    ### ENM_types, harmonic contains [O-O] already
    for (I, prot) in enumerate(Proteins)
        for Domain in UnfoldedRegions[prot]
            for i in range(Domain[1], Domain[2]-1)
                if !( (i,i+1) in BackboneDiffLengths[prot])
                    push!(BB_groups, (offsets[I]+i-1,offsets[I]+i)) ## shift from 1 -> 0, because of julia -> python
                end
            end
        end
    end
    return (N_Bonds+ENM_Bonds, ENM_types, vcat(BB_ID, ENM_typeid), vcat(BB_groups, ENM_groups), harmonic )
end

@doc raw"""
    ComputeHOOMD_ENMIndices(ConstraintDict, BackboneCorrectionDict, Sequences, Proteins)

Calculate the Indices, that are nessesary to creat a start file im HOOMD, .
    
**Arguments**
- `ConstraintDict::Dictionary with atoms, bond length and id which will be connected via ENM.
- `BackboneCorrectionDict`::Dictionary with atoms, bond length and id which will be connected via backbone.
- `Sequences`: The Sequences of the Proteins.
- `Proteins`: List of Protein Names.

**Return**:
- Number of bonds
- Vector of bond type names 
- ID connecting bond tuple to bond type
- Vector of tuples defining all bonds
- Dict{String, Dict{Symbol, Float64}} which defines the bonds as used in HOOMD.
"""
function ComputeHOOMD_ENMIndices(ConstraintDict, BackboneCorrectionDict , Sequences, Proteins)
    offsets = cumsum(length.(Sequences))

    bondLength = 0.38 # in nm
    BB_N = 0
    ENM_N = 0
    BB_ID = Int[]
    ENM_ID = Int[]
    BB_groups = Vector{Tuple{Int32, Int32}}()
    ENM_groups = Vector{Tuple{Int32, Int32}}()
    harmonic = Dict{String, Dict{Symbol, Float64}}()
    harmonic["O-O"] = Dict(:r => bondLength, :k => 8033.0) # default backbone bond values
    off = 0

    for (i, Prot) in enumerate(Proteins)
        for (atom_1, atom_2, r0, ind) in BackboneCorrectionDict[Prot]
            if r0 != bondLength
                push!(BB_ID, ind) 
            elseif r0 == bondLength
                push!(BB_ID, 0)
            end
            push!(BB_groups, (atom_1 -1 + off, atom_2 -1 + off)) ## shift from 1 -> 0, because of julia -> python
            BB_N+= 1
        end
        off = offsets[i]
    end
    off = 0
    for (i, Prot) in enumerate(Proteins)
        for (atom_1, atom_2, r0, ind) in ConstraintDict[Prot]
            push!(ENM_ID, ind) 
            push!(ENM_groups, (atom_1 -1 + off, atom_2 -1 + off))
            ENM_N+= 1
        end
        off = offsets[i]
    end
    BB_tmp  = vcat([BackboneCorrectionDict[prot] for prot in Set(Proteins)]...)
    ENM_tmp = vcat([ConstraintDict[prot]           for prot in Set(Proteins)]...)

    BB_types  = ["BB_$(ind)"  for ind in sort(collect(Set(getindex.(BB_tmp ,4)))) if ind != 0]
    ENM_types = ["ENM_$(ind)" for ind in sort(collect(Set(getindex.(ENM_tmp,4))))]

    for (r,ind) in zip(getindex.(BB_tmp,3) , getindex.(BB_tmp,4))
        if ind!=0 ### "BB_0" is the same as "O-O"
            harmonic["BB_$(ind)" ] = Dict(:r => r, :k => 8033)
        end
    end
    for (r,ind) in zip(getindex.(ENM_tmp,3), getindex.(ENM_tmp,4))
        harmonic["ENM_$(ind)" ] = Dict(:r => r, :k => 700)
    end
    
    return (BB_N + ENM_N, vcat(["O-O"], BB_types, ENM_types), vcat(BB_ID, ENM_ID), vcat(BB_groups, ENM_groups), harmonic)
end

@doc raw"""
    DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0, pae_cut=1.85)
Return a dictionary of atoms and there distances that are nessesary for the Elastic Network Model, and a dictionary for the correction of the backbone bonds in the ENM domain.
    
**Arguments**
- `BasePath`: The base path of the simulation setup procedure.
- `DomainDict`:: The Domains in which the ENM is active.
- `Proteins`: List of Protein Names.
- `ProteinJSON`: Dictionary of AlphaFold data of the Proteins in JSON format.
- `BBProtein`: The atom from which the AlphaFold datas are set for the aminoacid.
- `rcut`: Cut of length for the ENM in Angstroem.
- `plDDTcut`: Cut of plDDT parameter of AlphaFold reference for the ENM.
- `pae_cut`: Cut of pae parameter of AlphaFold reference for the ENM.

**Return**:
* `ConstraintDict`: Dictionary that maps Protein names to a Vector containing Tuples of Indices of i,j and distance r which define additional bonds necessary for ENM.
* `BackboneCorrectionDict`: Dictionary that maps Protein names to a Vector containing Tuples of Indices of i,j and distance r which define the backbone bonds that may have different lengths in ENM regions.
"""
function DetermineCalvados3ENMfromAlphaFold(BasePath::String, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0, pae_cut=1.85, precision=1000)
    ## distances in nm
    ciffolder = "$(BasePath)/InitFiles/CifFiles"
    ConstraintDict = Dict{String, Vector{Tuple{Int,Int, Float64, Int32}}}() ### Contains the bonds from ENM
    BackboneCorrectionDict = Dict{String, Vector{Tuple{Int,Int, Float64, Int32}}}() ### Contains the bonds from backbone with different lengths in folded regions
    bondLength = 0.38  ## coordinates are in nm

    ExistingConstraints=Dict{Float64, Int32}(bondLength=>1) ### maps existing constraint determined by distance onto integer
    ExistingBackbone   =Dict{Float64, Int32}(bondLength=>1) ### maps existing backbone   determined by distance onto integer
    cnt=1
    for Prot in Set(Proteins)
        if length(DomainDict[Prot])>0
            CifPath = "$(ciffolder)/$(Prot).cif"
            lines = readlines(CifPath)
            x = zeros(length(lines))
            y = zeros(length(lines))
            z = zeros(length(lines))
            plDDT = zeros(length(lines))

            step = 1
            for line in readlines(CifPath)
                fields = strip.(split(line))
                if !isempty(fields) && fields[1] == "ATOM" && fields[4] == BBProtein
                    x[step] = parse(Float64, fields[11])
                    y[step] = parse(Float64, fields[12])
                    z[step] = parse(Float64, fields[13])
                    plDDT[step]= parse(Float64, fields[15])
                    step += 1
                end
            end
            ConstraintDict[Prot] = []
            BackboneCorrectionDict[Prot] = []
            pae =JSON.parsefile(ProteinJSON[Prot])["pae"]
            for i in 1:step-1 # up to NAtom-1
                #in_any_domain = false
                for Domain in DomainDict[Prot]
                    if Domain[1] ≤ i ≤ Domain[2] && i + 1 ≤ Domain[2]
                        dist = sqrt((x[i+1]-x[i])^2 + (y[i+1]-y[i])^2 + (z[i+1]-z[i])^2)
                        dist = round(dist*precision)/precision ### reduces distance in bond length to .3f precision
                        dist /= 10.0 ### convert to nm
                        if !(dist in keys(ExistingBackbone))
                            ExistingBackbone[dist]=cnt
                            cnt+=1
                        end
                        push!(BackboneCorrectionDict[Prot], (i, i+1,dist, ExistingBackbone[dist]))

                        #in_any_domain = true
                        if plDDT[i] ≥ plDDTcut
                            for j in i+3:Domain[2]
                                dist_sqr = (x[j]-x[i])^2 + (y[j]-y[i])^2 + (z[j]-z[i])^2
                                if dist_sqr < rcut^2 && plDDT[j] ≥ plDDTcut && pae[i][j] < pae_cut && pae[j][i] < pae_cut
                                    dist=sqrt(dist_sqr)
                                    dist = round(dist*precision)/precision ### reduces distance in bond length to .3f precision
                                    dist /= 10.0 ### convert to nm
                                    if !(dist in keys(ExistingConstraints))
                                        ExistingConstraints[dist]=cnt
                                        cnt+=1
                                    end
                                    push!(ConstraintDict[Prot], (i, j, dist, ExistingConstraints[dist]))
                                end
                            end
                        end
                        break
                    end
                end
            
                #if !in_any_domain
                #    push!(BackboneCorrectionDict[Prot], (i, i+1, bondLength))
                #end
            end
        end
    end
    return ConstraintDict, BackboneCorrectionDict
end

@doc raw"""
    GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)

Generates a dictionary that defines the unfolded regions of proteins based on the definition of folded regions.
    
**Arguments**
- `Proteins`: List of Protein Names.
- `DomainDict`:: The Domains in which the ENM is active.
- `Sequences`: The Sequences of the Proteins.

**Return**:
- UnfoldedDict:: Dict{String, Vector{Tuple{Int64, Int64}}}(): Maps Protein names onto vectors of tuples that define the regions of unfolded regions.
"""
function GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)
    ProtLength = Dict([prot=>length(seq) for (seq, prot) in zip(Sequences, Proteins)])
    UnfoldedDict= Dict{String, Vector{Tuple{Int64, Int64}}}()

    for Prot in Set(Proteins)
        if length(DomainDict[Prot])>0
            N = ProtLength[Prot]

            FoldedDomains = sort(DomainDict[Prot])
            UnfoldedDomains = []
            if FoldedDomains[1][1]!=1
                push!(UnfoldedDomains, (1, FoldedDomains[1][1]))
            end
            for (i, val) in enumerate(FoldedDomains[1:end-1])
                push!(UnfoldedDomains, (val[2], FoldedDomains[i+1][1]))
            end
            if FoldedDomains[end][2]!=N
                push!(UnfoldedDomains, (FoldedDomains[end][2], N))
            end
            UnfoldedDict[Prot] = UnfoldedDomains
        end
    end
    return UnfoldedDict
end

@doc raw"""
    writeGSDStartFile(FileName::String, NAtoms::I, NBonds::I, NAngles::I, NDihedrals::I,Box::Vector{R}, Positions::Array{R}, AaToId::Dict{Char, <:Integer},Sequences,  InputImage::Array{I2}, InputMasses::Array{<:Real}, InputCharges::Array{R}, DihedralMap::Dict, DihedralList::Matrix{<:Integer}, AaToSigma::Dict{Char, <:Real}, UseAngles::Bool, SimulationType::String, Domains, ENM) where {I<:Integer, R<:Real, I2<:Integer}

A GSD data file is written, that include the parameters for the Simulation witch are given as Arguments.
    
**Arguments**
- `FileName::String`: Name of the output GSD file.
- `NAtoms::Int`: The total number of amino acids (atoms) in the system.
- `NBonds::Int`: The number of bonds between amino acids (`NAtoms - 1`).
- `NAngles::Int`: The number of angles formed between amino acids (`NAtoms - 2`).
- `NDihedrals::Int`: The number of dihedral angles between amino acids (`NAtoms - 3`).
- `Box::Vector{Float}`: 3-element vector specifying the dimensions of the simulation box
- `Positions::Array{Float}`: Array of the atomic coordinates of the Aminoacids.
- `AaToId::Dict{Char, Int}`: Dictionary mapping each amino acid type to a unique ID number.
- `Sequences::Array{String}`: A list of protein sequences, where each sequence is represented as a string of amino acids.
- `InputImage::Array{Int}`: An array used to determine periodic boundary conditions and correct atomic positions.
- `InputMasses::Array{Float}`: A 1D array specifying the mass of each amino acid.
- `InputCharges::Array{Float}`: A 1D array specifying the electric charge of each amino acid.
- `DihedralMap::Dict`: A dictionary mapping unique dihedral angle definitions (four atom indices) to dihedral types.
- `DihedralList::Matrix{Int}`: A matrix where each row defines a specific dihedral angle using atom indices.
- `AaToSigma::Dict{Char, <:Real}`: A dictionary mapping amino acid types to their Lennard-Jones σ-parameter (used in force field calculations).
- `UseAngles::Bool`: If `true`, angle and dihedral interactions are included in the GSD file.
- `SimulationType::String`: Type of simulation.
- `Domains::List(Int)`: Domains in which the ENM is active.
- `ENM::Tuple`: Data that are nessesary for ENM (bond name, id, group).

**Creates**:
* Writes the GSD data files.
"""
function writeGSDStartFile(FileName::String, NAtoms::I, NBonds::I, NAngles::I, NDihedrals::I,Box::Vector{R}, Positions::Array{R}, AaToId::Dict{Char, <:Integer},Sequences,  InputImage::Array{I2}, InputMasses::Array{<:Real}, InputCharges::Array{R}, DihedralMap::Dict, DihedralList::Matrix{<:Integer}, AaToSigma::Dict{Char, <:Real}, UseAngles::Bool, SimulationType::String, Domains, ENM) where {I<:Integer, R<:Real, I2<:Integer}
 
    snapshot = GSDFormat.Frame()    
    snapshot.configuration.step = 1 
    snapshot.configuration.dimensions = 3 
    snapshot.configuration.box = [Box[1],Box[2], Box[3], 0, 0, 0]./10.0 #4:6 are tilt
    snapshot.particles.N = NAtoms
    # pos = zeros(length(Sequences), maximum(length.(Sequences)), 3)
    # reshape -> [N_total, 3]
    snapshot.particles.position = reshape(permutedims(Positions, (2,1,3)), (size(Positions, 1)*size(Positions, 2), 3))./10.0 ### permute to get alignment in memory, reshape to match gsd formart
    IdToAa = Dict(value => key for (key, value) in AaToId)
    snapshot.particles.types =  [string(IdToAa[Id]) for Id in  sort(collect(values(AaToId))) ] #string.(collect(keys(AaToId)))
    snapshot.particles.typeid = [Int32(AaToId[AA])-1 for AA in join(Sequences)] ### convert to python numbering
    snapshot.particles.image = reshape(permutedims(InputImage, (2,1,3)), (size(InputImage, 1)*size(InputImage, 2), 3)) ### permute to get alignment in memory, reshape to match gsd formart
    snapshot.particles.mass = InputMasses
    snapshot.particles.charge = InputCharges
    snapshot.particles.diameter =  [Float32(AaToSigma[AA])/10.0  for AA in join(Sequences)]

    ### Bond_data.group = (self.N, getM(data)
    B_N, B_types, B_typeid, B_group_matrix = NBonds, ["O-O"], zeros(Int32, NBonds), getBonds(Sequences, M=2)
    # Create Bonds
    if SimulationType == "Calvados3"
        ENMB_N, ENMB_types, ENMB_typeid, ENMB_group_vector, harmonic = ENM
        B_N = ENMB_N
        B_types = ENMB_types
        B_typeid = ENMB_typeid
        ENMB_group_matrix = permutedims(hcat(collect.(ENMB_group_vector)...))
        B_group_matrix = ENMB_group_matrix
    end

    snapshot.bonds.N = B_N
    snapshot.bonds.types = B_types
    snapshot.bonds.typeid = B_typeid
    snapshot.bonds.group = B_group_matrix

    if UseAngles
        # Create Angles
        snapshot.angles.N = NAngles
        snapshot.angles.types = ["O-O-O"]
        snapshot.angles.typeid = zeros(Int32, NAngles)
        snapshot.angles.group = getBonds(Sequences, M=3)

        # Create Dihedrals
        snapshot.dihedrals.N =  NDihedrals 
        snapshot.dihedrals.types = ["$(ids[1])-$(ids[2])-$(ids[3])-$(ids[4])" for ids in collect(keys(DihedralMap))]#string.(collect(values(DihedralMap)))
        snapshot.dihedrals.typeid = [DihedralMap[DihedralList[key,:]]-1 for key in axes(DihedralList,1)] ### convert to python numbering
        snapshot.dihedrals.group = getBonds(Sequences, M=4)
    end

    file = GSDFormat.open(FileName, 'w')
    GSDFormat.append(file, snapshot)
    GSDFormat.close(file)
end

end
