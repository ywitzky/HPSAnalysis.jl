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
function computeInterChainContactMatrix(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0, NSubMatrices=5)::Matrix{R}  where {R<:Real, I<:Integer}
    ContactMatrices = [Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1])) for _ in 1:NSubMatrices]
    UsedChains = zeros(NSubMatrices)

    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    Sim.CellResolution = ones(3).*maximum(cutoffs)
    initCellLists(Sim)
    ### Elastic Network models screw up the estimation for the MaxParticlesPerCell
    Sim.CenterBox = zeros(R, (3 *Sim.MaxParticlesPerCell,3))
    Sim.NeighBox  = zeros(R, (27*Sim.MaxParticlesPerCell,3))
    Sim.CL_Dist   = zeros(R, (27*Sim.MaxParticlesPerCell^2,3)) ### 3 rows intermediate steps, first row contains the results

    ChainAtomsNeighbour = zeros(I, 27*Sim.MaxParticlesPerCell)
    ChainsNeighbour     = zeros(I, 27*Sim.MaxParticlesPerCell)

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)

    println("Range $(Sim.ClusterRange)")
    for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if step < Sim.EquilibrationTime || step %Sim.RGMeasureStep != 0
            continue 
        end ### ensures that there exist a clustering and while downsampling the evaluated frames
        println("Step $k of $(length(Sim.ClusterRange))")

        COM[Sim.SlabAxis] = AxisCOM[k]
        chains = isnothing(bounds) ? collect(1:Sim.NChains) : getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, R.(1.0./Sim.BoxLength))

        #UsedChains += length(chains)

        ### compute CellLists for step
        Sim.CellStep[1] = step
        resetCellLists(Sim)
        computeCellLists(Sim; ComputeCharges=false, ComputeChains=true) 

        ### iterate naively over all pairs of chains
        for (Cid, C) in enumerate(chains)
            C_subMatId = mod1(Cid, NSubMatrices)
            UsedChains[C_subMatId] += 1
            start = Sim.ChainStart[C]
            stop  = Sim.ChainStop[C]

            ### iterate over all cell lists where chain C is present
            for (xind, yind, zind) in getUniqueCellsOfChain(Sim, C)
                ChainAtomsInCenter = [id for id in Sim.CellList[xind, yind, zind] if id>=start && id <= stop]
                IMax = length(ChainAtomsInCenter)
                if IMax ==0 continue end
                ### get positions for center box
                Sim.CenterBox[1:IMax,1] .= Sim.x[ChainAtomsInCenter, step]
                Sim.CenterBox[1:IMax,2] .= Sim.y[ChainAtomsInCenter, step]
                Sim.CenterBox[1:IMax,3] .= Sim.z[ChainAtomsInCenter, step]

                ### have to compute all octants since we do not count all chains; aggregating all of them at once is quicker than individual computation -> SIMD 
                JMax = 1
                for xi in -1:1
                    xi_res = xind+xi
                    for yi in -1:1
                        yi_res = yind+yi
                        for zi in -1:1
                            zi_res = zind+zi
                            ### filter all particles that dont belong to chain C
                            for (i,id) in enumerate(Sim.CellList[xi_res, yi_res, zi_res])
                                if id<start || id > stop
                                    ChainAtomsNeighbour[JMax] = id 
                                    ChainsNeighbour[JMax] = Sim.ChainNumCellList[xi_res, yi_res, zi_res][i]
                                    JMax += 1
                                end
                            end
                        end
                    end
                end
                JMax -= 1
                            
                if JMax ==0 continue end

                ### get positions for comparison box
                Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour[1:JMax], step]
                Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour[1:JMax], step]
                Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour[1:JMax], step]

                ### compute distances in an optimized manor-> saved in Sim.CL_Dist[:,1]
                computeDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)

                cnt = 1
                for i in 1:IMax
                    atom_i = ChainAtomsInCenter[i]
                    rel_i = atom_i - start+1
                    cutoff_i = cutoffs[Sim.IDs[atom_i]]
                    #res_i = Sim.IDs[atom_i]
                    for j in 1:JMax
                        atom_j = ChainAtomsNeighbour[j] ### here we take neighbour cell

                        if Sim.CL_Dist[cnt,1]< (cutoff_i+ cutoffs[Sim.IDs[atom_j]])/2.0
                            rel_j = atom_j - Sim.ChainStart[ChainsNeighbour[j]] +1
                            ContactMatrices[C_subMatId][rel_i, rel_j] += 1
                        end
                        cnt += 1
                    end
                end

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


function computeInterChainMatrixHelper(Sim::SimData{R,I}, xind::I2, yind::I2, zind::I2, border::Bool, cutoffs::Vector{R2}, ValidChain::Vector{Bool}, maxcutff::R2) where {R<:Real,R2<:Real, I<:Integer,I2<:Integer}

    ChainAtomsInCenter =  [id for (C,id) in zip(Sim.ChainNumCellList[xind, yind, zind], Sim.CellList[xind,yind,zind]) if ValidChain[C]]
    ChainsInCenter     =  [C  for (C,id) in zip(Sim.ChainNumCellList[xind, yind, zind], Sim.CellList[xind,yind,zind]) if ValidChain[C]]

    IMax = length(ChainAtomsInCenter)
    if IMax ==0
        return
    end
    step = Sim.CellStep

    ### get positions for center box
    Sim.CenterBox[1:IMax,1] .= Sim.x[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]
    Sim.CenterBox[1:IMax,2] .= Sim.y[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]
    Sim.CenterBox[1:IMax,3] .= Sim.z[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]

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

                ChainAtomsNeighbour = Sim.CellList[xi_res,yi_res,zi_res]  #[id for (i,(C,id)) in enumerate(zip(Sim.ChainNumCellList[xi_res, yi_res, zi_res], Sim.CellList[xi_res,yi_res,zi_res])) ]
                JMax = length(ChainAtomsNeighbour)

                if JMax == 0 continue end
                ### get positions for center box
                Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour, step]
                Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour, step]
                Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour, step]

                ### compute distances in an optimized manor
                computeDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)

                for i in 1:IMax
                    for j in 1:JMax
                        if Sim.CL_Dist[cnt,1] <  maxcutff
                            atom_j = ChainAtomsNeighbour[j] ### here we take neighbour cell
                            C_j = Sim.ChainNumCellList[xi_res, yi_res, zi_res][j]
                            if ChainsInCenter[i] != C_j
                                atom_i = ChainAtomsInCenter[i]
                                rel_i = atom_i - Sim.ChainStart[ChainsInCenter[i]]+1
                                cutoff_i = cutoffs[Sim.IDs[atom_i]]
                                if Sim.CL_Dist[cnt,1]< (cutoff_i+ cutoffs[Sim.IDs[atom_j]])/2.0
                                    rel_j = atom_j - Sim.ChainStart[C_j]+1

                                    Sim.ContactMatrices[rel_i, rel_j] += 1
                                end
                            end
                        end
                        cnt += 1
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

    UsedChains = 0
    Sim.CellResolution = maximum(cutoffs)
    initCellLists(Sim)
    Sim.CenterBox = zeros(R, (Sim.MaxParticlesPerCell,3))
    Sim.NeighBox  = zeros(R, (Sim.MaxParticlesPerCell,3))
    Sim.CL_Dist   = zeros(R, (Sim.MaxParticlesPerCell^2,3)) ### 3 rows intermediate steps, first row contains the results

    x_low  = ceil(I,bounds[1][1]/Sim.CellResolution)
    x_high = ceil(I,bounds[1][2]/Sim.CellResolution) 
    y_low  = ceil(I,bounds[2][1]/Sim.CellResolution)
    y_high = ceil(I,bounds[2][2]/Sim.CellResolution) 
    z_low  = ceil(I,bounds[3][1]/Sim.CellResolution)
    z_high = ceil(I,bounds[3][2]/Sim.CellResolution) 

    ChainAtomsInCenter = Vector{I}(zeros(Sim.MaxParticlesPerCell))
    ChainsInCenter     = Vector{I}(zeros(Sim.MaxParticlesPerCell))

    chains = Vector{Bool}([true for _ in 1:Sim.NChains])

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
        computeInterChainMatrixHelper_(Sim::SimData{R,I}, xind::I2, yind::I2, zind::I2, border::Bool) where {R<:Real, I<:Integer,I2<:Integer}  = computeInterChainMatrixHelper(Sim, xind, yind, zind, border, cutoffs, chains, maximum(cutoffs))#, ChainAtomsInCenter, ChainsInCenter)

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