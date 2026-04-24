@doc raw"""
    Only created to debugging purposes use *computeInterChainContactMatrix* in production code. May be faster for very small systems.
"""
function computeInterChainContactMatrix_Naiv(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0) where {R<:Real, I<:Integer}
    #Sim.ContactMatrices = [zeros(R, Sim.ChainLength[C], Sim.ChainLength[C]) for C in 1:Sim.NChains]
    Sim.ContactMatrices = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))

    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    UsedChains = 0

    max_cutoff_sqr = maximum(cutoffs)^2

    for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        println("Step: $step")
        ### detect chains which are fully within bounds
        if isnothing(bounds)
            chains = collect(1:Sim.NChains)
        else
            chains = Vector{I}()
            for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
                (x_min, x_max) = extrema(Sim.x[start:stop,step])
                (y_min, y_max) = extrema(Sim.y[start:stop,step])
                (z_min, z_max) = extrema(Sim.z[start:stop,step])

                if  x_min>bounds[1][1] && x_max < bounds[1][2] && 
                    y_min>bounds[2][1] && y_max < bounds[2][2] && 
                    z_min>bounds[3][1] && z_max < bounds[3][2]
                    push!(chains,C)
                end
            end
        end
        UsedChains += length(chains)

        ### iterate naively over all pairs of chains
        for (Cid, C) in enumerate(chains)
            start = Sim.ChainStart[C]
            stop  = Sim.ChainStop[C]
            @inbounds for (i_rel,i) in enumerate(start:stop)
                i_cutoff = cutoffs[Sim.IDs[i]]

                ### periodic distance
                dist_x = abs.(Sim.x[:,step].-Sim.x[i,step]) 
                dist_y = abs.(Sim.y[:,step].-Sim.y[i,step]) 
                dist_z = abs.(Sim.z[:,step].-Sim.z[i,step]) 

                dist_x -= @. floor(I, dist_x * 1.0/Sim.BoxLength[1]+0.5)*Sim.BoxLength[1]
                dist_y -= @. floor(I, dist_y * 1.0/Sim.BoxLength[2]+0.5)*Sim.BoxLength[2]
                dist_z -= @. floor(I, dist_z * 1.0/Sim.BoxLength[3]+0.5)*Sim.BoxLength[3]

                dist_sqr = @. dist_x^2+dist_y^2+dist_z^2

                for j in findall(dist_sqr.<max_cutoff_sqr)
  
                    if  (j< start || j>stop) # C != C2
                        ChainOfj = searchsortedfirst(Sim.ChainStart, j) -1 
                        ChainOfj = ChainOfj != size(Sim.ChainStart,1) && Sim.ChainStart[ChainOfj+1]==j ? ChainOfj+1 : ChainOfj

                        j_rel = j - Sim.ChainStart[ChainOfj] +1

                        if dist_sqr[j]< ((cutoffs[Sim.IDs[j]]+i_cutoff)/2.0)^2
                            Sim.ContactMatrices[i_rel,j_rel] += 1
                        end
                    end
                end
            end
        end
    end
    println("UsedChains $UsedChains")
    Sim.ContactMatrices = (Sim.ContactMatrices.+transpose(Sim.ContactMatrices))./2.0
    Sim.ContactMatrices ./= UsedChains

    return Sim.ContactMatrices
end

@doc raw"""
    computeInterChainContactMatrix(Sim::SimData{R,I})

Computes the inter chain contact histrogram at each step of the *ClusterRange* by computing the distances between the residues of all proteins that are within the *bounds* in the *SlabAxis* and the neighbours the amino acids. The distance cutoffs for the histograms are defined by *rel_cutoff* ⋅ *CutoffDict*[type].

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `bounds::Union{Nothing, Vector{Tuple{R2, R2}}}`: Bound alongs the SlabAxis at which particles are considered to be in the bulk of the dense phase.
- `CutoffDict::Dict{Char,Float}=Sim.IDToSigmas`: Dictionary that defines the type specific factor for the cutoff.
- `rel_cutoff=1.0`: Global scaling factor for the cutoff.


**Create**:
* Sim.ContactMatrices
"""
function computeInterChainContactMatrix(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0, NSubMatrices=5,fac=5)::Matrix{R}  where {R<:Real, I<:Integer}
    ContactMatrices = [Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1])) for _ in 1:NSubMatrices]
    UsedChains = zeros(NSubMatrices)

    chain_of_atom = zeros(I, Sim.NAtoms)
    for C in 1:Sim.NChains
        for atom in Sim.ChainStart[C]:Sim.ChainStop[C]
            chain_of_atom[atom] = C
        end
    end


    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    Sim.CellResolution = ones(3).*maximum(cutoffs)
    initCellLists(Sim)
    ### Elastic Network models screw up the estimation for the MaxParticlesPerCell
    Sim.CenterBox = zeros(R, (fac*3 *Sim.MaxParticlesPerCell,3))
    Sim.NeighBox  = zeros(R, (fac*27*Sim.MaxParticlesPerCell,3))
    Sim.CL_Dist   = zeros(R, (fac*27*Sim.MaxParticlesPerCell^2,3)) ### 3 rows intermediate steps, first row contains the results

    ChainAtomsNeighbour = zeros(I, fac*27*Sim.MaxParticlesPerCell)
    ChainsNeighbour     = zeros(I, fac*27*Sim.MaxParticlesPerCell)

    ChainAtomsInCenter  = zeros(I, fac*5*27*Sim.MaxParticlesPerCell)

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)

    max_id = maximum(values(Sim.IDs))  # Assuming IDs are integers
    cutoff_combinations = zeros(R, max_id, max_id)
    for i in 1:max_id
        for j in 1:max_id
            cutoff_combinations[i,j] = ((cutoffs[i] + cutoffs[j])/2)^2
        end
    end

    invBoxLength= R.(1.0./Sim.BoxLength)
    for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if step < Sim.EquilibrationTime || step %Sim.RGMeasureStep != 0
            continue 
        end ### ensures that there exist a clustering and while downsampling the evaluated frames
        #println("Step $k of $(length(Sim.ClusterRange))")

        COM[Sim.SlabAxis] = AxisCOM[k]
        chains = isnothing(bounds) ? collect(1:Sim.NChains) : getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, invBoxLength)

        #UsedChains += length(chains)

        ### compute CellLists for step
        Sim.CellStep[1] = step
        resetCellLists(Sim)
        computeCellLists(Sim; ComputeCharges=false, ComputeChains=true) 

        #println("ASDFGH")
        ### iterate naively over all pairs of chains
        #chains=[471]
        for (Cid, C) in enumerate(chains)
          #  println("Cid $Cid, C $C")
            C_subMatId = mod1(Cid, NSubMatrices)
            UsedChains[C_subMatId] += 1
            start = Sim.ChainStart[C]
            stop  = Sim.ChainStop[C]

            ### iterate over all cell lists where chain C is present
           # println("A")
           # println(length(getUniqueCellsOfChain(Sim, C)))
            for (xind, yind, zind) in getUniqueCellsOfChain(Sim, C)
              #  println("xinds = $(xind) , $yind, $zind")
                IMax = 0
                @inbounds for id in Sim.CellList[xind, yind, zind]
                    if start <= id <= stop
                        IMax += 1
                        ChainAtomsInCenter[IMax] = id
                        Sim.CenterBox[IMax,1] = Sim.x[id, step]
                        Sim.CenterBox[IMax,2] = Sim.y[id, step]
                        Sim.CenterBox[IMax,3] = Sim.z[id, step]
                    end
                end
                #ChainAtomsInCenter = [id for id in Sim.CellList[xind, yind, zind] if id>=start && id <= stop]
                #IMax = length(ChainAtomsInCenter)
                if IMax ==0 continue end
                ### get positions for center box
                #Sim.CenterBox[1:IMax,1] .= Sim.x[ChainAtomsInCenter[1:IMax], step]
                #Sim.CenterBox[1:IMax,2] .= Sim.y[ChainAtomsInCenter[1:IMax], step]
                #Sim.CenterBox[1:IMax,3] .= Sim.z[ChainAtomsInCenter[1:IMax], step]


                ### have to compute all octants since we do not count all chains; aggregating all of them at once is quicker than individual computation -> SIMD 
                JMax = 0
                for zi in -1:1
                    zi_res = zind+zi
                    for yi in -1:1
                        yi_res = yind+yi
                        for xi in -1:1
                            xi_res = xind+xi
                            ### filter all particles that dont belong to chain C
                            #=
                            for (i,id) in enumerate(Sim.CellList[xi_res, yi_res, zi_res])
                                if id<start || id > stop
                                    ChainAtomsNeighbour[JMax] = id 
                                    ChainsNeighbour[JMax] = chain_of_atom[id]#Sim.ChainNumCellList[xi_res, yi_res, zi_res][i]
                                    JMax += 1
                                end
                            end
                            =#
                            cell   = Sim.CellList[xi_res, yi_res, zi_res]          # Vector{I}
                            chainv = Sim.ChainNumCellList[xi_res, yi_res, zi_res] # Vector{I} (same length)

                            # ---- fast loop over the cell (no enumerate, no bounds checks) -----
                            @simd for i in eachindex(cell)
                                id = cell[i]               # atom index
                                # “id not belonging to the current chain”
                                if id < start || id > stop
                                    JMax += 1
                                    ChainAtomsNeighbour[JMax] = id
                                    #ChainsNeighbour[JMax]    = chainv[id]   # or Sim.chain_of_atom[id] if you stored it globally
                                    
                                end
                            end
                        end
                    end
                end
                #println("B")

                if JMax ==0 continue end

                ### get positions for comparison box
                Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour[1:JMax], step]
                Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour[1:JMax], step]
                Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour[1:JMax], step]

                ### compute distances in an optimized manor-> saved in Sim.CL_Dist[:,1]
                #CB_tmp = @views Sim.CenterBox[1:IMax,:]
                #NB_tmp = @views Sim.NeighBox[1:JMax,:]
                computeSqrDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)

                cnt = 1
                for i in 1:IMax
                    #println("i $i")
                    atom_i = ChainAtomsInCenter[i]
                    rel_i = atom_i - start+1
                    cutoff_i = cutoffs[Sim.IDs[atom_i]]
                    id_i = Sim.IDs[ChainAtomsInCenter[i]]

                    #res_i = Sim.IDs[atom_i]
                    for j in 1:JMax
                        #println("j $j")
                        atom_j = ChainAtomsNeighbour[j] ### here we take neighbour cell
                        id_j = Sim.IDs[ChainAtomsNeighbour[j]]
                        #cutoff_combinations[id_i,id_j]  #
                        #println("$j, $id_i, $id_j, $cnt")
                        if Sim.CL_Dist[cnt,1]< cutoff_combinations[id_i,id_j]  #((cutoff_i+ cutoffs[Sim.IDs[atom_j]])/2.0)^2
                            #rel_j = atom_j - Sim.ChainStart[ChainsNeighbour[j]] +1
                            rel_j = atom_j - Sim.ChainStart[chain_of_atom[atom_j]] +1
                            ContactMatrices[C_subMatId][rel_i, rel_j] += 1
                        end
                        cnt += 1
                    end
                end
                #println("D")

            end
        end
    end

    mean = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))
    error = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))

    for i in 1:NSubMatrices
        println("UsedChains $(UsedChains[i])")
        ContactMatrices[i] = (ContactMatrices[i].+transpose(ContactMatrices[i]))./2.0
        ContactMatrices[i] ./= UsedChains[i]
        mean .+= ContactMatrices[i] 
    end
    mean ./= NSubMatrices

    for i in 1:NSubMatrices
        error .+= (mean .- ContactMatrices[i]).^2
    end
    error .= sqrt.(error)./sqrt(NSubMatrices)

    Sim.ContactMatrices = mean
    Sim.ContactMatricesError = error
end


function computeInterChainMatrixHelper(Sim::SimData{R,I}, xind::I2, yind::I2, zind::I2, border::Bool, cutoffs::Vector{R2}, SqrCutoffMatrix::Matrix{R}, ValidChain::Vector{Bool}, maxcutff::R2,chain_of_atom::Vector{I}, ChainAtomsInCenter::Vector{I},ChainsInCenter::Vector{I}, ChainAtomsNeighbour::Vector{I}) where {R<:Real,R2<:Real, I<:Integer,I2<:Integer}

    #ChainAtomsInCenter =  [id for (C,id) in zip(Sim.ChainNumCellList[xind, yind, zind], Sim.CellList[xind,yind,zind]) if ValidChain[C]]
    #ChainsInCenter     =  [C  for (C,id) in zip(Sim.ChainNumCellList[xind, yind, zind], Sim.CellList[xind,yind,zind]) if ValidChain[C]]
    step = Sim.CellStep[1]


    cell_chains = Sim.ChainNumCellList[xind, yind, zind]
    cell_atoms  = Sim.CellList[xind,yind,zind]
    IMax=0
    @inbounds @simd for i in eachindex(cell_atoms)   # avoids bounds checks & enables SIMD
        c = cell_chains[i]               # chain number of the i‑th atom in the cell
        if ValidChain[c]                # branch is SIMD‑friendly (mask)
            IMax += 1
            id = cell_atoms[i]
            ChainAtomsInCenter[IMax] = id    # atom id
            ChainsInCenter[IMax] = c         # chain id (you could also use `Sim.chain_of_atom[atom]`)
            Sim.CenterBox[IMax,1] = Sim.x[id, step]
            Sim.CenterBox[IMax,2] = Sim.y[id, step]
            Sim.CenterBox[IMax,3] = Sim.z[id, step]
        end
    end

    #IMax = length(ChainAtomsInCenter)
    if IMax ==0
        return
    end

    ### get positions for center box
    #Sim.CenterBox[1:IMax,1] .= Sim.x[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]
    #Sim.CenterBox[1:IMax,2] .= Sim.y[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]
    #Sim.CenterBox[1:IMax,3] .= Sim.z[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]

    JMax =0
    atom_i = 0
    atom_j = 0
    res_i=0
    #cnt = 1

    ### compute all octants since not every chain is taken
    for xi in -1:1
        xi_res = xind+xi
        for yi in -1:1
            yi_res = yind+yi
            for zi in -1:1
                zi_res = zind+zi
                cnt = 1

                #ChainAtomsNeighbour = Sim.CellList[xi_res,yi_res,zi_res]  #[id for (i,(C,id)) in enumerate(zip(Sim.ChainNumCellList[xi_res, yi_res, zi_res], Sim.CellList[xi_res,yi_res,zi_res])) ]
                #JMax = length(ChainAtomsNeighbour)

                #if JMax == 0 continue end
                ### get positions for center box
                #Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour, step]
                #Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour, step]
                #Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour, step]

                cell   = Sim.CellList[xi_res, yi_res, zi_res]          # Vector{I}
                chainv = Sim.ChainNumCellList[xi_res, yi_res, zi_res] # Vector{I} (same length)

                # ---- fast loop over the cell (no enumerate, no bounds checks) -----
                @simd for i in eachindex(cell)
                    id = cell[i]               # atom index

                        JMax += 1
                        ChainAtomsNeighbour[JMax] = id
                end
            end
        end
    end
    Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour[1:JMax], step]
    Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour[1:JMax], step]
    Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour[1:JMax], step]

    ### compute distances in an optimized manor
    computeSqrDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)

    cnt=0
    for i in 1:IMax
        @inbounds @simd  for j in 1:JMax
            cnt += 1
            if Sim.CL_Dist[cnt,1] <  maxcutff
                atom_j = ChainAtomsNeighbour[j] ### here we take neighbour cell
                C_j= chain_of_atom[atom_j]
                #C_j = Sim.ChainNumCellList[xi_res, yi_res, zi_res][j]
                if ChainsInCenter[i] != C_j
                    atom_i = ChainAtomsInCenter[i]
                    id_i= Sim.IDs[atom_i]
                    rel_i = atom_i - Sim.ChainStart[ChainsInCenter[i]]+1
                    cutoff_i = cutoffs[Sim.IDs[atom_i]]
                    id_j= Sim.IDs[atom_j]
                    if Sim.CL_Dist[cnt,1]< SqrCutoffMatrix[id_j,id_i] #(cutoff_i+ cutoffs[Sim.IDs[atom_j]])/2.0
                        rel_j = atom_j - Sim.ChainStart[C_j]+1

                        Sim.ContactMatrices[rel_i, rel_j] += 1
                    end
                end
            end
        end
    end

    return nothing
end


@doc raw"""
    Expected to be faster than *computeInterChainContactMatrix* but so far still slower for the tested system sizes.
"""
function computeInterChainContactMatrix_fast(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0::R2)::Matrix{R} where {R<:Real, I<:Integer} 
    #Sim.ContactMatrices = [zeros(R, Sim.ChainLength[C], Sim.ChainLength[C]) for C in 1:Sim.NChains]
    Sim.ContactMatrices = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))

    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    chain_of_atom = zeros(I, Sim.NAtoms)
    for C in 1:Sim.NChains
        for atom in Sim.ChainStart[C]:Sim.ChainStop[C]
            chain_of_atom[atom] = C
        end
    end

    UsedChains = 0
    Sim.CellResolution = ones(R,3).*maximum(cutoffs)
    initCellLists(Sim)
    Sim.CenterBox = zeros(R, (Sim.MaxParticlesPerCell,3))
    Sim.NeighBox  = zeros(R, (27*Sim.MaxParticlesPerCell,3))
    Sim.CL_Dist   = zeros(R, (27*Sim.MaxParticlesPerCell^2,3)) ### 3 rows intermediate steps, first row contains the results

    x_low  = ceil(I,bounds[1][1]/Sim.CellResolution[1])
    x_high = ceil(I,bounds[1][2]/Sim.CellResolution[1]) 
    y_low  = ceil(I,bounds[2][1]/Sim.CellResolution[2])
    y_high = ceil(I,bounds[2][2]/Sim.CellResolution[2]) 
    z_low  = ceil(I,bounds[3][1]/Sim.CellResolution[3])
    z_high = ceil(I,bounds[3][2]/Sim.CellResolution[3]) 

    ChainAtomsInCenter  = Vector{I}(zeros(I,Sim.MaxParticlesPerCell))
    ChainsInCenter      = Vector{I}(zeros(I,Sim.MaxParticlesPerCell))
    ChainAtomsNeighbour = Vector{I}(zeros(I,27*Sim.MaxParticlesPerCell))

    chains = Vector{Bool}([true for _ in 1:Sim.NChains])

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)

    #ChainAtomsInCenter  = zeros(I, 5*27*Sim.MaxParticlesPerCell)
    #ChainAtomsInCenter  = zeros(I, 5*27*Sim.MaxParticlesPerCell)

    max_id = maximum(values(Sim.IDs))  # Assuming IDs are integers
    cutoff_combinations = zeros(R, max_id, max_id)
    for i in 1:max_id
        for j in 1:max_id
            cutoff_combinations[i,j] = ((cutoffs[i] + cutoffs[j])/2)^2
        end
    end


    for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        println("Step $k of $(length(Sim.ClusterRange))")
        COM[Sim.SlabAxis] = AxisCOM[k]
        chains_tmp = getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, R.(1.0./Sim.BoxLength))
        chains.=true
        for id in chains_tmp
            chains[id]=true
        end
        UsedChains += sum(chains)

        ### compute CellLists for step
        Sim.CellStep[1] = step
        resetCellLists(Sim)
        computeCellLists(Sim; ComputeCharges=false, ComputeChains=true) 

        ### iterate over all pairs of chains
        computeInterChainMatrixHelper_(Sim::SimData{R,I}, xind::I2, yind::I2, zind::I2, border::Bool) where {R<:Real, I<:Integer,I2<:Integer}  = computeInterChainMatrixHelper(Sim, xind, yind, zind, border, cutoffs,cutoff_combinations, chains, maximum(cutoffs),chain_of_atom, ChainAtomsInCenter, ChainsInCenter,ChainAtomsNeighbour )#, ChainAtomsInCenter, ChainsInCenter)

        iterateThroughCellList(Sim, computeInterChainMatrixHelper_; CellMaxima=[x_low+1, x_high, y_low+1, y_high, z_low+1, z_high])
    end
    Sim.ContactMatrices = (Sim.ContactMatrices+transpose(Sim.ContactMatrices))./2.0
    Sim.ContactMatrices ./= UsedChains

    Sim.ContactMatrices
end

@doc raw"""
    computeIntraChainContactMatrix(Sim::SimData{R,I})

Computes the intra chain contact map by averaging over all pair distances first and taking the logarithm thereof. 

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* Sim.IntraChainContactMatrix
"""
function computeIntraChainContactMatrix(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.IntraChainContactMatrix = [zeros(R, Sim.ChainLength[C], Sim.ChainLength[C]) for C in 1:Sim.NChains]

    for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
            @inbounds for (i_rel,i) in enumerate(start:stop)
                for j in i+1:stop
                    j_rel = j-Sim.ChainStart[C]+1
                    dist_sqr = (Sim.x_uw[i,step]-Sim.x_uw[j,step])^2+(Sim.y_uw[i,step]-Sim.y_uw[j,step])^2+(Sim.z_uw[i,step]-Sim.z_uw[j,step])^2
                    dist = sqrt(dist_sqr)
                    Sim.IntraChainContactMatrix[C][i_rel,j_rel] += dist
                end
            end
        end
    end
    
    for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
        Sim.IntraChainContactMatrix[C] ./= length( Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        @inbounds for (i_rel,i) in enumerate(start:stop)
            for j in i+1:stop 
                j_rel = j-Sim.ChainStart[C]+1
                Sim.IntraChainContactMatrix[C][i_rel,j_rel] = log(Sim.IntraChainContactMatrix[C][i_rel,j_rel])
                Sim.IntraChainContactMatrix[C][j_rel,i_rel] = Sim.IntraChainContactMatrix[C][i_rel,j_rel]
            end
        end
    end

    Sim.IntraChainContactMatrix
end