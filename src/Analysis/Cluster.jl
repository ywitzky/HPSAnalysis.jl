export computeClusterCOMs

@doc raw"""
    computeCOMClusters(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=50.0) where{T<:Real, Int<:Integer}

Computes a List of Cluster based on the vicinity of the COMs.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `Cutoff::Float`: How close the chains have to be for a Cluster.

**Returns**:
- `Cluster::Vector{Vector{Vector{Int}}}`: List of Clusters of the chains.
"""
function computeCOMClusters(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=50.0) where{T<:Real, Int<:Integer}
    ### cutoff from dignon et al. Sequence determinants of protein phase behavior from a coarse-grained model
    Clusters = Vector{Vector{Vector{Int}}}()

    pre = Sim.ChargeAnalysisCutoff 
    Sim.ChargeAnalysisCutoff = Cutoff
    initCellLists(Sim)
    for (i,step) in enumerate(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        G = Graphs.SimpleGraph{Int}()
        add_vertices!(G, Sim.NChains)

        for C in 1:Sim.NChains
            dx = Sim.COM[C+1:Sim.NChains,1, step] .- Sim.COM[C,1, step]
            dy = Sim.COM[C+1:Sim.NChains,2, step] .- Sim.COM[C,2, step]
            dz = Sim.COM[C+1:Sim.NChains,3, step] .- Sim.COM[C,3, step]

            dist_sqr = dx.^2 .+dy.^2 .+dz.^2
            for (k,K) in enumerate(C+1:Sim.NChains)
                if dist_sqr[k]<Cutoff^2
                    add_edge!(G, C,K)
                end
            end
        end

        fill!(Sim.ChargeChainContactMatrix, 0)
        push!(Clusters,  weakly_connected_components(G))

    end
    Sim.ChargeAnalysisCutoff = pre


    Sim.ClusterRange = Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
    Sim.Clusters = Clusters
    return Clusters
end

function computeChargeClusters(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=5.0) where{T<:Real, Int<:Integer}
    printstyled("Function computeClusters does only detect charge contacts for clusters.\n"; color=:yellow)
    ### abuse cellList mechanism developed for Charge contacts
    Sim.ChargeResidueContactMatrix = zeros(eltype(Sim.x), Sim.MaxChainLength,Sim.MaxChainLength)
    Sim.ChargeChainContactMatrix = zeros(eltype(Sim.x), Sim.NChains,Sim.NChains)
    Clusters = Vector{Vector{Vector{Int}}}()

    pre = Sim.ChargeAnalysisCutoff 
    Sim.ChargeAnalysisCutoff = Cutoff
    initCellLists(Sim)
    for (i,step) in enumerate(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        Sim.ChargeContactValid = Dict()
        
        G = Graphs.SimpleGraph{Int}()
        add_vertices!(G, Sim.NChains)
        
        Sim.CellStep[1] = step
        HPSAnalysis.resetCellLists(Sim) ### empty all cells
        if i%100==0
            println("Step: $(Sim.CellStep) of $(Sim.NSteps)")
        end
        HPSAnalysis.computeCellLists(Sim) ### compute CellLists  for step

        iterateThroughCellList(Sim, computeChargeBondTimesForCell )

        for key in keys(Sim.ChargeContactValid)
            atom_i = key[1] ### always positive charge
            atom_j = key[2] ### always negative charge

            ### compute charge bond Matrix

            ### compute chain numbers
            chain_num_i = Int32(ceil(atom_i/Sim.MaxChainLength))
            chain_num_j = Int32(ceil(atom_j/Sim.MaxChainLength))

            ### compute residue number in chain
            id_in_chain_i = atom_i - (chain_num_i-1)*Sim.MaxChainLength
            id_in_chain_j = atom_j - (chain_num_j-1)*Sim.MaxChainLength
            
            ### increase Contactmatrices
            Sim.ChargeResidueContactMatrix[id_in_chain_i, id_in_chain_j] += 1
            Sim.ChargeChainContactMatrix[chain_num_i, chain_num_j] += 1
        end

        for C in 1:Sim.NChains
            for K in C+1:Sim.NChains
                if Sim.ChargeChainContactMatrix[C,K]>0
                    add_edge!(G, C,K)
                end
            end
        end

        fill!(Sim.ChargeChainContactMatrix, 0)
        push!(Clusters,  weakly_connected_components(G))

    end
    Sim.ChargeAnalysisCutoff = pre


    Sim.ClusterRange = Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
    Sim.Clusters = Clusters
    return Clusters
end

function computeClusterCOMs(Sim::SimData{T,I} ) where {I<:Integer, T<:Real}
    Coms= Vector{Matrix{T}}()
    for (Cstep,C) in enumerate(Sim.Clusters)
        step = Sim.ClusterRange[Cstep]
        tmp = zeros(length(C),3 )
        for (CID,Cluster) in enumerate(C)
            mass=0
            for id in Cluster
                tmp[CID,:].+= Sim.COM[id,:,step]*Sim.ChainMasses[CID]
                mass+= Sim.ChainMasses[CID]
            end
            tmp[CID,:]./=mass
        end
        push!(Coms, tmp)
    end
    return Coms
end
