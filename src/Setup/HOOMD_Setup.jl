function WriteHOOMDSequences(filename::String, Sequences)
    io = open(filename, "w");

    for (SeqID, Seq) in enumerate(Sequences)
        write(io, "$Seq\n")
    end
    close(io);
end

function WriteDictionaries(filename::String, ToCharge, ToID, ToMass, ToDiameter, ToLambda)
    io = open(filename, "w")
    write(io, "// ID,resname, Charge, Mass, λ   \n")

    for key in sort(collect(keys(ToID)), by=x->ToID[x])
        write(io, " $(ToID[key]), $key, $(ToCharge[key]), $(ToMass[key]), $(ToDiameter[key]), $(ToLambda[key])\n")
    end
    close(io);
end

function WriteHOOMDParticlesInput(filename::String, pos::Array{R}, ToCharge, ToID, Sequences, ToMass, ToDiameter, image) where {R<:Real}

    io = open(filename, "w");
    write(io, "### N, id, x ,y ,z , charge, m, diameter \n");
    cnt=1
    for (SeqID, Seq) in enumerate(Sequences)
        for (atom,res) in enumerate(Seq)
            write(io, "$(cnt), $(ToID[res] -1) , $(@sprintf("%.3f",pos[SeqID,atom, 1])), $(@sprintf("%.3f",pos[SeqID,atom, 2])), $(@sprintf("%.3f",pos[SeqID,atom, 3])), $(ToCharge[res]), $(ToMass[res]), $(ToDiameter[res]) , $(image[SeqID, atom,1]) , $(image[SeqID, atom,2]) , $(image[SeqID, atom,3]) \n");
            cnt += 1
        end
    end

    close(io);
end

function WriteDihedrals(filename, dihedral_map, dihedral_eps)
    io = open(filename, "w");
    write(io , "# $(dihedral_eps) \n")
    for key in keys(dihedral_map)
        write(io, "$(key[1]) $(key[2]) $(key[3]) $(key[4]) $(dihedral_map[key])\n")
    end

    close(io)
end

function WriteParams(filename, SimName, Temp, NSteps, NOut, Timestep, Box, Seed; Minimise=true, TrajectoryName="traj.gsd", UseAngles=true, UseCharge=true, Alt_GSD_Start="-", Create_Start_Config=false, ϵ_r=1.73136, κ=1.0, Device="GPU")
    io = open(filename, "w");
    write(io, "Simname: $SimName\n")
    write(io, "Seed: $Seed\n")
    write(io, "Temp: $Temp\n")
    write(io, "NSteps: $NSteps\n")
    write(io, "NOut: $NOut\n")
    write(io, "dt: $Timestep\n")
    write(io, "Lx: $(Box[1])\n")
    write(io, "Ly: $(Box[2])\n")
    write(io, "Lz: $(Box[3])\n")
    write(io, "Minimise: $(Minimise)\n")
    write(io, "Trajectory: $(TrajectoryName)\n")
    write(io, "UseAngles: $(UseAngles)\n")
    write(io, "UseCharge: $(UseCharge)\n")
    write(io, "Alt_GSD_Start: $(Alt_GSD_Start)\n")
    write(io, "Create_Start_Config: $(Create_Start_Config)\n")
    write(io, "epsilon_r: $(ϵ_r)\n")
    write(io, "kappa: $(κ)\n")
    write(io, "Device: $(Device)\n")

    close(io);
end
