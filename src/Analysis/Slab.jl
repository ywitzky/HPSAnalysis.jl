@doc raw"""
    computeSlabHistogram(Sim::SimData{R,I}; Use_Alpha=false, Use_Types=false) where {R<:Real, I<:Integer}
Computes centered slab histograms along Sim.SlabAxis.

Centered slab histograms are computed for all amino acids, only the positive and only the negatives at all times. **Use\_Alpha** enables the computation of alpha helices and their own slab histogram. Amino acids specific histograms are enabled through **Use\_Types**. 

Results are not return but stored in Sim.SlabHistogramSeries as an Offset array where the first dimension ranges from -boxwidth/Sim.Resolution:Sim.Resolution:boxwidth/Sim.Resolution. Default Sim.Resolution is set to ``1~\AA``. The second index are the steps at which clusters and slab histogram were computed according to Sim.ClusterRange. The third index are the amino acids which have been used the order in which they are mentioned above. If used, amino acid specific histograms have the index 4+Sim.IDs. 

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Creat**:
* `Sim.SlabHistogramSeries`: Stores mass densities across slabs for different atom types.
"""
function computeSlabHistogram(Sim::SimData{R,I}; Use_Alpha=false, Use_Types=false) where {R<:Real, I<:Integer}
    if Use_Alpha && sum(Sim.TorsionAngles[:,1])==0 
        computeDihedralAngles(Sim) ### ensure that TorsionAngles have been computed
    end

    AxisCOM = computeCOMOfLargestCluster(Sim)

    if Sim.SlabAxis==1
        SlabCoord = Sim.x_uw
    elseif Sim.SlabAxis==2
        SlabCoord = Sim.y_uw
    elseif Sim.SlabAxis == 3
        SlabCoord = Sim.z_uw
    else 
        ArgumentError("SlabAxis is not properly specified.")
    end

    Len = deepcopy(Sim.BoxLength[Sim.SlabAxis])
    Len_inv = 1/Len

    NHists = 1 + 2 +1 + Sim.NAtomTypes*Use_Types # one normal, one for positive charge, one for negative charge, one for alpha helices and one for each type

    ### array with 1:N steps for -boxwidth:boxwitdh in the direction of the slab
    Sim.SlabHistogramSeries = OffsetArray(zeros(R, Int32(ceil((Sim.BoxSize[Sim.SlabAxis,2]-Sim.BoxSize[Sim.SlabAxis,1])/Sim.Resolution)) , Int32(Sim.NSteps), NHists), Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) , 1:Sim.NSteps, 1:NHists)


    AllAxis = [1,2,3]
    deleteat!(AllAxis, Sim.SlabAxis)


    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    volume =((Sim.BoxSize[AllAxis[1],2]-Sim.BoxSize[AllAxis[1],1])*(Sim.BoxSize[AllAxis[2],2]-Sim.BoxSize[AllAxis[2],1]))
    conversion = 1.66053906660/volume/Sim.Resolution

    if Use_Alpha
        Pseudohelical = zeros(Bool,Sim.NAtoms)
        AlphaHelixProb = zeros(eltype(Sim.x), Sim.NAtoms)
    end

    for (j,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if j %200 == 0 println("step $j/$(length(Sim.ClusterRange))") end

        if Use_Alpha ### compute helical states
            Pseudohelical .= false
            AlphaHelixProb .= 0
            decidePseudoHelicals(Sim, Pseudohelical, AlphaHelixProb, step)
        end

        for atom in 1:Sim.NAtoms 
            pos = (SlabCoord[atom,step]-AxisCOM[j]-0.5) ### center slab at 0
            pos -= Len*round(I, pos*Len_inv)            ### wrap back to central images
            ind = ceil(Int32,((pos)/Sim.Resolution))    ### get index according to resolution

            Sim.SlabHistogramSeries[ind , step,1]+= Sim.Masses[atom]
            if Sim.Charges[atom] > 0
                Sim.SlabHistogramSeries[ind , step,2]+= Sim.Masses[atom]
            elseif Sim.Charges[atom] < 0
                Sim.SlabHistogramSeries[ind , step,3]+= Sim.Masses[atom]
            end

            if Use_Alpha && AlphaHelixProb[atom]>0
                Sim.SlabHistogramSeries[ind, step, 4] += Sim.Masses[atom]
            end

            if Use_Types
                Sim.SlabHistogramSeries[ind, step, Sim.IDs[atom]+4]+= Sim.Masses[atom] ### lowest ID is 1
            end
        end
    end
    Sim.SlabHistogramSeries *= conversion 

    nothing
end

function computeDensityHistogram(Sim::HPSAnalysis.SimData{R,I}, DivLength=I(10))  where {R<:Real, I<:Integer}
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

@doc raw"""
    computeSlabHistogram(Sim::SimData; Width=25, MaxVal=0.9, Surface_fac=0.8)

Computes average density within dense phase and dilute phase as well as the indices below/above which Sim.SlabHistogramSeries is in dense/dilute phase.

A mirror symmetric density around zero is computed from which dense phase approximation ρ\_app is defined as the mean value of the first **Width** steps. The dense phase boundary is **Surface_fac** times the distance **r_dense** at which the density drops below **MaxVal** times ρ\_app. Similarly the dilute phase boundary is the distance at which the density drops below **1-MaxVal** times ρ\_app plus **r_dense** times (1-**Surface_fac**).
"""
function computeSlabDensities(Sim::HPSAnalysis.SimData{R,I}; Width=25, MaxVal=0.9, Surface_fac=0.8)  where {R<:Real, I<:Integer}

    ### compute mirror symmetric density of centered histogram
    N = extrema(axes(Sim.SlabHistogramSeries,1))[2]-1
    density = zeros(R,N+1)
    NSteps = size(Sim.SlabHistogramSeries,2)

    ρ_dense  = zeros(R,NSteps)
    ρ_dilute = zeros(R, NSteps)
    ind_dense  = zeros(I, NSteps)
    ind_dilute = zeros(I, NSteps)

    for step in axes(Sim.SlabHistogramSeries,2)
        fill!(density,zero(R))

        for i in 0:N 
            density[i+1] += Sim.SlabHistogramSeries[i, step,1]
            density[i+1] += Sim.SlabHistogramSeries[-i, step,1]
        end
        density /= 2.0

        avg = mean(density[1:Width])

        tmp = findlast(density.>avg*MaxVal)
        diff = tmp *(1-Surface_fac)
        ind_dense[step] = round(I, tmp*Surface_fac)
        ρ_dense[step] = mean(density[1:ind_dense[step]])

        N_dilute = 2*ind_dense[step]*1/Surface_fac< N ? 2*ind_dense[step]*1/Surface_fac :  N-Width
        N_dilute =  ceil(I, N_dilute)
        avg = mean(density[N_dilute:end])

        ind_dilute[step] = ceil(I,findfirst(density.-avg .<(1.0-MaxVal)*(ρ_dense[step]-avg)) + diff )

        ρ_dilute[step] = mean(density[ind_dilute[step]:end])

        ### ### convert to indexing starting at 0 for usage in offset array of SlabHistogramSeries
        ind_dense[step] -= 1
        ind_dilute[step] -= 1 
    end

    return ρ_dense, ρ_dilute, ind_dense, ind_dilute
end