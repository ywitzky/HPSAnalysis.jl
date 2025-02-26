# Setup for a Simulation
First, the setup begins by creating a data structur that organizes the information and Simulation datas, for each protein and temperture. For the proteins there corresponding sequence is taken from a dictionary.\
Based on the protein sequence, the number of chains (NChains) required to achieve the target density in the simulation box is calculated. Once NChains is determined, a list is generated containing NChains repetitions of the protein sequences. \
Next, the parameters of the simulation box, such as its size and boundary conditions, are defined. \
Using the previously defined data, the initial coordinates of the proteins are calculated, ensuring a proper starting configuration for the simulation. \
Finally, all necessary inputs are written, including the simulation parameters, the interaction model (Yukawa potential with Debye-Hückel screening), relevant dictionaries, and the start file. These files will be used to initialize the simulation.


```julia
ToCreate =  ["RS31", "RS31a_IDR2", "RS31_IDR2"]

for (protID, protein) in enumerate(ToCreate)
    mkpath(BasePath*"$(protein)/")
    for temp in Temperatures
        ###Generation of the folder structur for each protein and temperature
        pad = lpad(getkey(RunsPerProtein,protein,1),3,"0")
        mkpath(BasePath*"$(protein)/$(temp)K/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/Restart/")
        Path = BasePath*"$(protein)/$(temp)K/RUN_$(pad)/"
        cd(Path)

        ###Sequence from predefined sequences 
        Seq = HPS.ProteinSequences.NameToSeq[protein]

        ###Calculation of number of chains (NChains) required to achieve the target density (Box=[L,L*factor,L])
        NChains = DetermineNumberChains(Seq, SimulationType="Calvados2", pH=pH, SideLength=SideLength, density=0.4,fac=width_multiplier)

        ###List of NChains repetitions of sequences and proteins
        Seq = HPS.ProteinSequences.NameToSeq[protein]
        Sequences= [deepcopy(Seq) for _ in 1:NChains]
        Proteins = [deepcopy(protein) for _ in 1:NChains]
        Info ="SLAB Simulation script for $protein.\n\n"

        ###Parameters of the simulation box
        BoxLengthShort=Float32(350.0)      
        BoxLengthLong=Float32(1500.)
        BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]

        SimName = "$(protein)_$temp"

        ###Initial coordinates of the proteins are calculated
        (pos, Data) = HPS.CreateStartConfiguration(SimName,Path , Float32.([BoxLengthShort,BoxLengthShort*width_multiplier , BoxLengthShort]), Proteins, Sequences, Regenerate=false; Axis="y")#Erstellung der Start Conformation

        ###Inputs are written, including the simulation parameters, the interaction model, relevant dictionaries, and the start file
        HPS.Setup.writeStartConfiguration("./$(protein)_slab","./$(SimName)_Start_slab.txt", Info, Sequences, BoxSize , 300_000_000, HOOMD=true, ; SimulationType="Calvados2" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH)
    end
end
```