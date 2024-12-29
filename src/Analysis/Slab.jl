function computeSlabHistogram(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    ### ensure that TorsionAngles have been computed
    if sum(Sim.TorsionAngles[:,1])==0
        computeDihedralAngles(Sim)
    end

    Clust_Coms = computeClusterCOMs(Sim)
    LClustID = argmax.([length.(C) for C in Sim.Clusters]) ### returns index for clusters with most proteins inside

    AxisCOM = [Clust_Coms[i][id, Sim.SlabAxis] for (i,id) in enumerate(LClustID)]


    #println(length(AxisCOM))
    if Sim.SlabAxis==1
        SlabCoord = Sim.x
    elseif Sim.SlabAxis==2
        SlabCoord = Sim.y  
    elseif Sim.SlabAxis == 3
        SlabCoord = Sim.z 
    else 
        ArgumentError("SlabAxis is not properly specified.")
    end


    NHists = 1 + 2 +1 + Sim.NAtomTypes # one normal, one for positive chage, one for negative charge, one for alpha helices and one for each type

    ### array with 1:N steps for -boxwidth:boxwitdh in the direction of the slab
    Sim.SlabHistogramSeries = OffsetArray(zeros(eltype(Sim.x), Int32(ceil((Sim.BoxSize[Sim.SlabAxis,2]-Sim.BoxSize[Sim.SlabAxis,1])/Sim.Resolution)) , Int32(Sim.NSteps), NHists), Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) , 1:Sim.NSteps, 1:NHists)

    AllAxis = [1,2,3]
    deleteat!(AllAxis, Sim.SlabAxis)


    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    volume =((Sim.BoxSize[AllAxis[1],2]-Sim.BoxSize[AllAxis[1],1])*(Sim.BoxSize[AllAxis[2],2]-Sim.BoxSize[AllAxis[2],1]))
    conversion = 1.66053906660/volume/Sim.Resolution
    lowestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1
    highestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) 
    Ninds = abs(lowestind)+1+highestind
    #printstyled("Subtract COM and add wrap afterwards\n"; color=:red)

    Pseudohelical = zeros(Bool,Sim.NAtoms)
    AlphaHelixProb = zeros(eltype(Sim.x), Sim.NAtoms)

    println(Sim.BoxSize)
    println("low: $(lowestind), high: $(highestind)")
    for (j,step) in enumerate(Sim.ClusterRange)### 1:NSteps
        if step %200 == 0 println("step $step") end

        ### compute helical states
        Pseudohelical .= false
        AlphaHelixProb .= 0
        decidePseudoHelicals(Sim, Pseudohelical, AlphaHelixProb, step)

        for atom in 1:Sim.NAtoms ### 1:Atoms
            ### so far dont use particles that arent in the box
            ind = Int32(ceil((SlabCoord[atom,step]-AxisCOM[j])/Sim.Resolution))


            #ind = ind >= lowestind && ind<= highestind ?  ind : continue
            if ind < lowestind 
                ind = highestind - (lowestind-ind)%Ninds
            elseif  ind> highestind
                ind =lowestind + (ind-highestind)%Ninds
            end

            #= 
            if ind < lowestind || ind> highestind
                println("atom $(atom)")
                println("COM: $(AxisCOM[j])")
                println("y: $(SlabCoord[atom,step])")
                println("ind pre: $(Int32(ceil((SlabCoord[atom,step]-AxisCOM[j])/Sim.Resolution)))")
                println("ind: $(ind)")
                println("bla: $(ind%Ninds)")
                println("fix: $( highestind - (lowestind-(ind%Ninds)))")
                println("fix2: $( highestind - (lowestind-ind)%Ninds)")
            end=#



            Sim.SlabHistogramSeries[ind , step,1]+= Sim.Masses[atom]
            if Sim.Charges[atom] > 0
                Sim.SlabHistogramSeries[ind , step,2]+= Sim.Masses[atom]
            elseif Sim.Charges[atom] < 0
                Sim.SlabHistogramSeries[ind , step,3]+= Sim.Masses[atom]
            end

            
            if AlphaHelixProb[atom]>0
                Sim.SlabHistogramSeries[ind, step, 4] += Sim.Masses[atom]
            end


            Sim.SlabHistogramSeries[ind, step, Sim.IDs[atom]+4]+= Sim.Masses[atom] ### lowest ID is 1

        end
        for ind in axes(Sim.SlabHistogramSeries,1)
            for type in axes(Sim.SlabHistogramSeries,3)
                Sim.SlabHistogramSeries[ind, step,type] *= conversion
            end
        end 
    end
end

function centerSlabHistogram(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    ### kann weg eigentlich....

    SlabCoord  = Sim.COM[:,Sim.SlabAxis,:]
    COM_inds = zeros(eltype(Sim.NSteps), Sim.NSteps)
    TotalMass= sum(Sim.ChainMasses)

    for step in 1:Sim.NSteps
        COM_inds[step] = Int32(ceil((sum(SlabCoord[:,step].*Sim.ChainMasses)/TotalMass)/Sim.Resolution))
    end
    lowestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1
    highestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) 

    for (step, ind) in enumerate(COM_inds)
        if ind > 0 
            tmp = Sim.SlabHistogramSeries[lowestind:lowestind+ind,step,:]
            Sim.SlabHistogramSeries[lowestind:highestind-ind, step,:] .= Sim.SlabHistogramSeries[lowestind+ind:highestind, step,:]
            Sim.SlabHistogramSeries[highestind-ind:highestind, step,:] .= tmp
        elseif ind < 0
            ind  = abs(ind)
            tmp = Sim.SlabHistogramSeries[highestind-ind:highestind,step,:]
            Sim.SlabHistogramSeries[lowestind+ind:highestind, step,:] .= Sim.SlabHistogramSeries[lowestind:highestind-ind, step,:]
            Sim.SlabHistogramSeries[lowestind:lowestind+ind, step,:] .= tmp
        end
    end
end

function computeDensityHistogram(Sim::LammpsAnalysis.SimData{R,I}, DivLength=I(10))  where {R<:Real, I<:Integer}
    if !(Sim.SlabAxis in [1,2,3])
        ArgumentError("SlabAxis is not properly specified.")
    end

    dims = DivLength*ones(I, 3)
    dims[Sim.SlabAxis] *= DivLength### subdivide long axis for finer grid
    offset  = Sim.BoxSize[:,1]
    divider = Sim.BoxLength./R.(dims)
    BoxHist = zeros(R,  (I.(dims)...) ) ### atom , step in x,y,z makes this more efficient since memory is aligned in last index
    indices = zeros(I,  (I.(dims)...) ) 
    ind = zeros(I, 3)
    tmp = zeros(R, 3)

    bools = zeros(Bool,3)
    bools_tmp = zeros(Bool,3)

    cnt = 0 

    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    volume = prod(divider)
    NRes = 1000
    conversion = R(1.66053906660/volume*NRes) ### *NRes for indicing later
    Sim.DensityHist = zeros(R, 2*NRes)
    xoff = I(dims[3]*dims[2])
    yoff = I(dims[3])
    for step in  Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        for atom in 1:Sim.NAtoms
            
            ind[1] = ceil(I,(Sim.x[atom,step]-offset[1])/ divider[1])
            ind[2] = ceil(I,(Sim.y[atom,step]-offset[2])/ divider[2])
            ind[3] = ceil(I,(Sim.z[atom,step]-offset[3])/ divider[3])
            
            bools[1] = ind[1]<=0 || ind[1]>dims[1]
            bools[2] = ind[2]<=0 || ind[2]>dims[2]
            bools[3] = ind[3]<=0 || ind[3]>dims[3]
            
            if bools[1] || bools[2] || bools[3]
                cnt += 1
                continue
            end

            @inbounds BoxHist[ind[1], ind[2], ind[3]] += Sim.Masses[atom]
        end
        BoxHist .*=  conversion
        indices .= ceil.(I,BoxHist).+I(1)

        for i in indices
            Sim.DensityHist[i] += R(1)
        end
    end

    Sim.DensityHist ./= length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)*prod(dims)
end