using Base.Iterators, HDF5

function getSlabCoordinate(Sim::SimData{R,I};Unwrapped=true) where {R<:Real, I<:Integer}
    if Sim.SlabAxis==1
        SlabCoord = Unwrapped ? Sim.x_uw : Sim.x
    elseif Sim.SlabAxis==2
        SlabCoord = Unwrapped ? Sim.y_uw : Sim.y
    elseif Sim.SlabAxis == 3
        SlabCoord = Unwrapped ? Sim.z_uw : Sim.z
    else 
        ArgumentError("SlabAxis is not properly specified.")
    end
    return SlabCoord
end

@inline function getRecenteredPositions(coord, com , Len, Len_inv)
    pos = (coord-com)                       ### center slab at 0
    pos -= Len*round(Int32, pos*Len_inv)    ### wrap back to central images
end

@inline function getRecenteredPositions(SlabCoord::Array{R}, atom ,j, step, AxisCOM, Len, Len_inv) where {R<:Real}
    getRecenteredPositions(SlabCoord[atom,step],AxisCOM[j], Len, Len_inv )
end

@inline function getChainsInBounds(Sim::SimData{R,I}, bounds::Union{Nothing, Vector{Tuple{R, R}}}, COM::Vector{R}, step, Len::Vector{R}, Len_inv::Vector{R}) where {R<:Real, I<:Integer}
    chains = Vector{I}()
    if !isnothing(bounds)
        for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
            (x_min, x_max) = getRecenteredPositions.(extrema(Sim.x_uw[start:stop, step]), COM[1], Len[1], Len_inv[1])
            (y_min, y_max) = getRecenteredPositions.(extrema(Sim.y_uw[start:stop, step]), COM[2], Len[2], Len_inv[2])
            (z_min, z_max) = getRecenteredPositions.(extrema(Sim.z_uw[start:stop, step]), COM[3], Len[3], Len_inv[3])

            if  x_min>bounds[1][1] && x_max < bounds[1][2] && 
                y_min>bounds[2][1] && y_max < bounds[2][2] && 
                z_min>bounds[3][1] && z_max < bounds[3][2] 

                push!(chains,C)
            #else 
            #    println("$y_min, $y_max, $(bounds[2])")
            end
        end
    end
    return chains
end

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
function computeSlabHistogram(Sim::SimData{R,I}; Ranges::Union{Nothing, Vector{UnitRange{I}}}=nothing, Use_Alpha=false, Use_Types=false) where {R<:Real, I<:Integer}
    if Use_Alpha && sum(Sim.TorsionAngles[:,1])==0 
        computeDihedralAngles(Sim) ### ensure that TorsionAngles have been computed
    end

    AxisCOM = computeCOMOfLargestCluster(Sim)

    SlabCoord =  getSlabCoordinate(Sim)

    Len = deepcopy(Sim.BoxLength[Sim.SlabAxis])
    Len_inv = 1.0/Len

    # all residue types, one for positive charged, one for negative charged, one for alpha helices,  one for each type, and according to range definition
    NHists = 1 + 2 +1 + Sim.NAtomTypes*Use_Types 
    rangeoffset = 1 + 2 +1 + Sim.NAtomTypes*Use_Types 
    if !isnothing(Ranges)
        NHists += length(Ranges)
    end

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

    lowestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1
    highestind =Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution))

    for (j,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if j %200 == 0 println("step $j/$(length(Sim.ClusterRange))") end

        if Use_Alpha ### compute helical states
            Pseudohelical .= false
            AlphaHelixProb .= 0
            decidePseudoHelicals(Sim, Pseudohelical, AlphaHelixProb, step)
        end

        for C in 1:Sim.NChains
            for (atom_in_chain, atom) in enumerate(Sim.ChainStart[C]:Sim.ChainStop[C])
                pos = getRecenteredPositions(SlabCoord, atom,j , step, AxisCOM, Len, Len_inv)
                ind = ceil(Int32,((pos-0.5)/Sim.Resolution))    ### get index according to resolution
                if ind == lowestind-1
                    ind = highestind
                end

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

                if !isnothing(Ranges)
                    atom_in_chain .∈ Ranges 
                    for (k, in_range) in enumerate(atom_in_chain .∈ Ranges)
                        if in_range
                            Sim.SlabHistogramSeries[ind, step, rangeoffset+k]+= Sim.Masses[atom]
                        end
                    end
                end
            end
        end
    end
    Sim.SlabHistogramSeries *= conversion 

    nothing
end

@doc raw"""
    computeDensityHistogram(Sim::SimData{R,I}, DivLength=I(10)) where {R<:Real, I<:Integer}
Computes a logarithmic histogram of densities of the subcubes of the simulation box.

The simulation box is divided into *DivLength*^4 subboxes where each dimension is divided into *DivLength* many subsections and the axes according to Sim.SlabAxis is divided into *DivLength*^2 many subsections.

Subboxes that are empty or those that contain the boundary of the simulation box are not considered for the histogram.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `DivLength=I(10))`: Division factor determining the length of each subbox. 

**Creat**:
* `Sim.DensityHist`: Stores subcube densities with logarithmic indices.
"""
function computeDensityHistogram(Sim::HPSAnalysis.SimData{R,I}, DivLength=I(10))  where {R<:Real, I<:Integer}
    if !(Sim.SlabAxis in [1,2,3])
        ArgumentError("SlabAxis is not properly specified.")
    end

    dims = DivLength*ones(I, 3)
    dims[Sim.SlabAxis] *= DivLength### subdivide long axis for finer grid
    offset  = Sim.BoxSize[:,1]
    divider = Sim.BoxLength./R.(dims)
    BoxHist = zeros(R,  (I.(dims)...) ) ### atom , step in x,y,z makes this 
    ind = zeros(I, 3)

    bools = zeros(Bool,3)
    cnt = 0 

    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    volume = prod(divider)
    res = 80 ### take 80 minorticks per power of 10
    conversion = R(1.66053906660/volume)
    Sim.DensityHist = zeros(R, 8*res) ### overexpect indices from 10^-6 to 10^1
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
        indices = ceil.(I,log10.(filter(x->x!=0.0 ,BoxHist[:])).*res).+I(6*res)

        for i in indices
            Sim.DensityHist[i] += R(1)
        end
        fill!(BoxHist, 0.0)
    end
    Sim.DensityHist ./= length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)*prod(dims)
end

function computeSlabDensities_Help!(SlabHistogramSeries::OffsetArrays.OffsetArray{R, 3, Array{R, 3}}, density::Array{R},step::I) where {R<:Real, I<:Int}
    fill!(density,zero(R))
    N = extrema(axes(SlabHistogramSeries,1))[2]-1

    for i in 0:N 
        density[i+1] += SlabHistogramSeries[i, step,1]
        density[i+1] += SlabHistogramSeries[-i, step,1]
    end
    density ./= 2.0
    return density
end



@doc raw"""
    computeSlabHistogram(Sim::SimData; Width=25, MaxVal=0.9, Surface_fac=0.8)

Computes average density within dense phase and dilute phase as well as the indices below/above which Sim.SlabHistogramSeries is in dense/dilute phase.

A mirror symmetric density around zero is computed from which dense phase approximation ρ\_app is defined as the mean value of the first **Width** steps. The dense phase boundary is **Surface_fac** times the distance **r_dense** at which the density drops below **MaxVal** times ρ\_app. Similarly the dilute phase boundary is the distance at which the density drops below **1-MaxVal** times ρ\_app plus **r_dense** times (1-**Surface_fac**).
"""
function computeSlabDensities(Sim::HPSAnalysis.SimData{R,I}; Width=25, MaxVal=0.9, Surface_fac=0.8, safety=75.0, global_inds=(nothing, nothing))  where {R<:Real, I<:Integer}

    safety_ = ceil(Int32, safety/Sim.Resolution)

    ### compute mirror symmetric density of centered histogram
    N = extrema(axes(Sim.SlabHistogramSeries,1))[2]-1
    density = zeros(R,N+1)
    NSteps = size(Sim.SlabHistogramSeries,2)

    ρ_dense  = zeros(R,NSteps)
    ρ_dilute = zeros(R, NSteps)
    ind_dense  = zeros(I, NSteps)
    ind_dilute = zeros(I, NSteps)

    #if all(isnothing.(global_inds)) ### default automatic setting of dense/dilute phases
        for step in axes(Sim.SlabHistogramSeries,2)
            computeSlabDensities_Help!(Sim.SlabHistogramSeries, density, step)

            avg = mean(density[1:Width])

            tmp = findlast(density.>avg*MaxVal)
            diff = tmp *(1-Surface_fac)
            ind_dense[step] = round(I, tmp*Surface_fac)
            ρ_dense[step] = mean(density[1:ind_dense[step]]) ### will be overwritten later

            N_dilute = 2*ind_dense[step]*1/Surface_fac< N ? 2*ind_dense[step]*1/Surface_fac :  N-Width
            N_dilute =  ceil(I, N_dilute)
            avg = mean(density[N_dilute:end])

            ind_dilute[step] = ceil(I,findfirst(density.-avg .<(1.0-MaxVal)*(ρ_dense[step]-avg)))
        end

        
        global_ind_dilute = isnothing(global_inds[1]) ? min(maximum(ind_dilute)+safety_, N) : global_inds[1]
        global_ind_dense  = isnothing(global_inds[2]) ? minimum(ind_dense) : global_inds[2]
        
    #else  ### manually set phase borders
    #    global_ind_dilute = global_inds[1]
    #    global_ind_dense  = global_inds[2]
    #end

    for step in axes(Sim.SlabHistogramSeries,2)
        computeSlabDensities_Help!(Sim.SlabHistogramSeries, density, step)
        ρ_dense[step] = mean(density[1:global_ind_dense]) 
        ρ_dilute[step] = mean(density[global_ind_dilute:end])

        ### ### convert to indexing starting at 0 for usage in offset array of SlabHistogramSeries
        #ind_dense[step] -= 1
        #ind_dilute[step] -= 1 

        ind_dense[step]  = global_ind_dense - 1
        ind_dilute[step] = global_ind_dilute - 1 
    end

    return ρ_dense, ρ_dilute, ind_dense, ind_dilute
end

@inline function getVoxelIndex(pos, res, off)
    tmp = ceil(Int32,((pos)/res))+off
    #tmp > res ? tmp-2*off :  (tmp< 1 ? tmp+2*off : tmp)
end


function computeBinderCumulantsOfSlabDensities(Sim::HPSAnalysis.SimData{R,I}, indices_dilute, indices_dense)  where {R<:Real, I<:Integer}

    #DenseCumulantBoxes, DiluteCumulantBoxes = computeBinderCumulantsSubBoxes(Sim, indices_dilute, indices_dense)

    order_param = Sim.SlabDenseCumulantBoxes .- Sim.SlabDiluteCumulantBoxes

    Sim.BinderSlabCumulants = zeros(5)
    norm = zeros(5)

    for (wid,width) in enumerate([1,2,3,4,6])
        for ax1 in partition(1:12, width)
            for ax2 in partition(1:12, width)
                tmp = order_param[ax1,ax2,:].^2
                Sim.BinderSlabCumulants[wid] += mean(tmp)^2/mean(tmp.^2)
                norm[wid] += 1
            end
        end
    end
    Sim.BinderSlabCumulants ./= norm
end

function computeBinderCumulantsSubBoxes(Sim::HPSAnalysis.SimData{R,I}, indices_dilute, indices_dense)  where {R<:Real, I<:Integer}
    AxisCOM = computeCOMOfLargestCluster(Sim)

    AllAxis = [1,2,3]
    deleteat!(AllAxis, Sim.SlabAxis)

    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    d1 = Sim.BoxLength[AllAxis[1]]/12.0
    d2 = Sim.BoxLength[AllAxis[2]]/12.0
    volume =d1*d2
    conversion = 1.66053906660/volume/Sim.Resolution

    ### get Long Axis and short axis depending on selection
    SlabCoord =  getSlabCoordinate(Sim;Unwrapped=false )
    tmp = deepcopy(Sim.SlabAxis)
    Sim.SlabAxis= AllAxis[1]
    Axis1 =  getSlabCoordinate(Sim;Unwrapped=false)
    Sim.SlabAxis= AllAxis[2]
    Axis2 =  getSlabCoordinate(Sim;Unwrapped=false)
    Sim.SlabAxis= tmp

    Len = deepcopy(Sim.BoxLength[Sim.SlabAxis])
    Len_inv = 1.0/Len

    ### divide short axes into 12 sub boxes
    DenseCumulantBoxes  = zeros(R, 12,12, Sim.NSteps)
    DiluteCumulantBoxes = zeros(R, 12,12, Sim.NSteps)

    dilute_cutoff = indices_dilute .* Sim.Resolution
    dense_cutoff  = indices_dense  .* Sim.Resolution
    for (j,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if j %200 == 0 println("step $j/$(length(Sim.ClusterRange))") end

        for atom in 1:Sim.NAtoms
            ### recenter histogram according to Cluster COMs
            pos = getRecenteredPositions(SlabCoord, atom ,j, step, AxisCOM, Len, Len_inv)

            if abs(pos)< dense_cutoff[step]
                ind1 = getVoxelIndex(Axis1[atom,step], d1, 6)
                ind2 = getVoxelIndex(Axis2[atom,step], d2, 6)
                DenseCumulantBoxes[ind1, ind2, j] += Sim.Masses[atom]
            end

            if abs(pos)> dilute_cutoff[step]
                ind1 = getVoxelIndex(Axis1[atom,step], d1, 6)
                ind2 = getVoxelIndex(Axis2[atom,step], d2, 6)
                DiluteCumulantBoxes[ind1, ind2, j] += Sim.Masses[atom]
            end
        end
        DenseCumulantBoxes[:,:,j] /= dense_cutoff[step]*2 
        DiluteCumulantBoxes[:,:,j] /= (Sim.BoxLength[Sim.SlabAxis]-dilute_cutoff[step]*2 )
    end
    Sim.SlabDenseCumulantBoxes  = DenseCumulantBoxes .*conversion
    Sim.SlabDiluteCumulantBoxes = DiluteCumulantBoxes.*conversion

    return Sim.SlabDenseCumulantBoxes, Sim.SlabDiluteCumulantBoxes
end

@doc raw"""
    computeSurfaceTension(Sim::HPSAnalysis.SimData{R,I},start::I;filename::String = "pressure.h5",tau::I = I(1),NSub::I= I(10))

Computes the **mechanical (pressure‑tensor) surface tension** of a slab
simulation that has been written out by HOOMD‑Blue.  The routine reads the
pressure‑tensor time series from an HDF5 file, extracts the normal and
tangential components with respect to the slab normal (`Sim.SlabAxis`),
applies the standard formula  

\[
\gamma = \frac{L_{\mathrm{slab}}}{2}\,
\bigl\langle P_{norm.} - P_{tang.} \bigr\rangle ,
\]

and returns the mean value together with an estimate of its statistical
uncertainty obtained by block‑averaging.

**Arguments**

- `Sim::HPSAnalysis.SimData{R,I}`:  A simulation data structure containing the Simulation information.
- `start::I`: index (in multiples of `Sim.TrajWriteOutFreq`) of the  first frame to be included in the analysis.

**Keyword arguments**

- `filename::String = "pressure.h5"`: name of the HDF5 file that contains the pressure tensor (`hoomd-data/md/compute/ThermodynamicQuantities/pressure_tensor`).
- `tau::I = I(1)`: stride for sub‑sampling the pressure‑tensor series
  (e.g. `tau=10` uses every 10‑th stored frame).
- `NSub::I = I(10)`: number of sub‑samples used for the block‑averaging
  error estimate.

**Returns**

A tuple `(γ, Δγ)` where  

- `γ` (`Sim.SurfaceTension`) is the mean surface tension,
- `Δγ` (`Sim.SurfaceTensionError`) is the estimated statistical error."""
function computeSurfaceTension(Sim::HPSAnalysis.SimData{R,I},start::I;filename::String="pressure.h5", tau::I=I(1), NSub::I=I(10)) where {R<:Real, I<:Integer}
    longfilename = Sim.BasePath*filename 

    file = h5open(longfilename) 
    pressure_tensor = file["hoomd-data/md/compute/ThermodynamicQuantities/pressure_tensor"]
    timestep = file["hoomd-data/Simulation/timestep"]

    ### start is assumed to be in multiples of the frequency of the frame outputs
    ### pressure data is written out way more

    start_time  = start * Sim.TrajWriteOutFreq
    start_index = findfirst(x->(x)>=start_time, timestep[:] )
    
    #https://hoomd-blue.readthedocs.io/en/v6.0.0/hoomd/md/compute/thermodynamicquantities.html#hoomd.md.compute.ThermodynamicQuantities.pressure_tensor
    p_xx = pressure_tensor[1, start_index:tau:end]
    p_yy = pressure_tensor[4, start_index:tau:end]
    p_zz = pressure_tensor[6, start_index:tau:end]

    all_axis = [p_xx, p_yy, p_zz]
    normal_p = all_axis[Sim.SlabAxis]
    tangential_p = sum(deleteat!(all_axis,Sim.SlabAxis))./2.0
    
    γ = Sim.BoxLength[Sim.SlabAxis]/2.0 .*(normal_p.-tangential_p)
    Sim.SurfaceTension = mean(γ)
    n_measure = floor(I, length(normal_p)/NSub)
    y_mean_sub = [mean(sub) for sub in Iterators.partition(γ, n_measure)][1:end-1] ### last partition is not full...
    Sim.SurfaceTensionError = sqrt(sum((y_mean_sub.-Sim.SurfaceTension).^2)/NSub)

    return Sim.SurfaceTension, Sim.SurfaceTensionError
end