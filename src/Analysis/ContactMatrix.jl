
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
                Sim.IntraChainContactMatrix[C][j_rel,i_rel]  = Sim.IntraChainContactMatrix[C][i_rel,j_rel]
            end
        end
    end
    

    Sim.IntraChainContactMatrix
end