function initCellLists(Sim::LammpsAnalysis.SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.CellStep = zeros(I, 1)
    Sim.CellDimensions= zeros(3,2)
    for i in 1:3
        for j in 1:2
            Sim.CellDimensions[i,j] = Int32(ceil(Sim.BoxSize[i,j]/Sim.CellResolution))
        end
    end

    CellVolume=Sim.CellResolution^3
    Sim.MaxParticlesPerCell= Int32(ceil(CellVolume/(4/3*pi*(4.5/2.)^3) )+1)# volume of the smallest residues

    xdim = 3+Sim.CellDimensions[1,2]-Sim.CellDimensions[1,1]
    ydim = 3+Sim.CellDimensions[2,2]-Sim.CellDimensions[2,1]
    zdim = 3+Sim.CellDimensions[3,2]-Sim.CellDimensions[3,1]
    xrange = (Sim.CellDimensions[1,1]-1):(Sim.CellDimensions[1,2]+1) 
    yrange = (Sim.CellDimensions[2,1]-1):(Sim.CellDimensions[2,2]+1)
    zrange = (Sim.CellDimensions[3,1]-1):(Sim.CellDimensions[3,2]+1)

    Sim.CellList = OffsetArray(Array{Vector{eltype(Sim.NSteps)}}(undef, xdim, ydim, zdim), xrange , yrange, zrange)
    Sim.PositiveCellList = OffsetArray(Array{Vector{eltype(Sim.NSteps)}}(undef, xdim, ydim, zdim), xrange , yrange, zrange)
    Sim.NegativeCellList = OffsetArray(Array{Vector{eltype(Sim.NSteps)}}(undef, xdim, ydim, zdim), xrange , yrange, zrange)

    resetCellLists(Sim)
end

function resetCellLists(Sim::SimData{T,I}) where {T<:Real, I<:Integer}  
    for zi in axes(Sim.CellList,3)
        for yi in axes(Sim.CellList,2)
            for xi in axes(Sim.CellList,1)
                Sim.CellList[xi, yi, zi] =  Vector{I}()
                Sim.PositiveCellList[xi, yi, zi] =  Vector{I}()
                Sim.NegativeCellList[xi, yi, zi] =  Vector{I}()
            end
        end
    end
end

function destroyCellLists(Sim::SimData{T,I}) where {T<:Real, I<:Integer}
    Sim.CellLists= zeros(0)
    Sim.PositiveCellLists= zeros(0)
    Sim.NegativeCellLists= zeros(0)
end

function computeCellLists(Sim::SimData{T,I} ) where {T<:Real, I<:Integer}
    step = Sim.CellStep[1]
    resetCellLists(Sim)

    x_low  = ceil(I,Sim.BoxSize[1,1]/Sim.CellResolution)
    x_high = ceil(I,Sim.BoxSize[1,2]/Sim.CellResolution) 
    y_low  = ceil(I,Sim.BoxSize[2,1]/Sim.CellResolution)
    y_high = ceil(I,Sim.BoxSize[2,2]/Sim.CellResolution) 
    z_low  = ceil(I,Sim.BoxSize[3,1]/Sim.CellResolution)
    z_high = ceil(I,Sim.BoxSize[3,2]/Sim.CellResolution) 

    xdiff =convert(I,1+x_high-x_low)
    ydiff =convert(I,1+y_high-y_low)
    zdiff =convert(I,1+z_high-z_low)

    ind  = 1
    xind = 1
    yind = 1
    zind = 1
    
    @views xind_arr = ceil.(I,Sim.x[:,step]/Sim.CellResolution)
    @views yind_arr = ceil.(I,Sim.y[:,step]/Sim.CellResolution)
    @views zind_arr = ceil.(I,Sim.z[:,step]/Sim.CellResolution)
    
    atom_cnt = 0
    for atom in axes(Sim.x,1) #1:Sim.NAtoms
            
        xind = xind_arr[atom]
        xind = xind >= x_low ? xind :  xind+xdiff
        xind = xind <= x_high ? xind : xind-xdiff

        yind = yind_arr[atom]
        yind = yind >= y_low  ? yind : yind+ydiff
        yind = yind <= y_high ? yind : yind-ydiff
        
        zind = zind_arr[atom]
        zind = zind >= z_low  ? zind : zind+zdiff
        zind = zind <= z_high ? zind : zind-zdiff

        ### all atoms in cell 
        push!(Sim.CellList[xind, yind, zind], atom)

        ### all negative residues
        if Sim.Charges[atom]<0
            push!(Sim.NegativeCellList[ xind, yind, zind], atom)
        end

        ### all positives residues
        if Sim.Charges[atom]>0
            push!(Sim.PositiveCellList[xind, yind, zind], atom)
        end
    end

    ### copy borders on the other side...

    ### iterate over celllists
    for  CellList in [Sim.CellList,Sim.NegativeCellList,Sim.PositiveCellList]

        for yind in y_low-1:y_high+1
            for zind in z_low-1:z_high+1
                for e in CellList[x_low, yind, zind]
                    ### copy -x to x and back
                    if ~(e in CellList[x_high+1, yind, zind])
                        push!(CellList[x_high+1, yind, zind], e )
                    end

                    ### copy x to -x and back
                    if ~(e in CellList[ x_low-1, yind, zind])
                        push!( CellList[x_low-1, yind, zind], e)
                    end
                end
            end
        end

        for xind in x_low-1:x_high+1
            for zind in z_low-1:z_high+1
                for e in CellList[xind, y_low, zind]
                    ### copy -y to y and back
                    if ~(e in CellList[xind , y_high+1, zind])
                        push!( CellList[ xind, y_high+1, zind], e)
                    end

                    ### copy y to -y and back
                    if ~(e in CellList[xind, y_low-1, zind])
                        push!(CellList[ xind, y_low-1, zind], e)
                    end
                end
            end
        end


        for xind in x_low-1:x_high+1
            for yind in y_low-1:y_high+1
                for e in CellList[xind, yind, z_low]
                    ### copy -z to z and back
                    if ~(e in CellList[xind , yind, z_high+1])
                        push!(CellList[xind, yind, z_high+1], e)
                    end

                    ### copy z to -z and back
                    if ~(e in CellList[xind, yind, z_low-1])
                        push!(CellList[xind, yind, z_low-1], e)
                    end
                end
            end
        end
    end
end

function iterateThroughCellList(Sim::SimData{T,I}, func::Function)where {T<:Real, I<:Integer}
    step = Sim.CellStep[1]
    x_low  = convert(Int,Sim.CellDimensions[1,1])
    x_high = convert(Int,Sim.CellDimensions[1,2])
    y_low  = convert(Int,Sim.CellDimensions[2,1])
    y_high = convert(Int,Sim.CellDimensions[2,2])
    z_low  = convert(Int,Sim.CellDimensions[3,1])
    z_high = convert(Int,Sim.CellDimensions[3,2])

    if x_high-x_low < 3 || y_high-y_low < 3 || z_high-z_low < 3
        printstyled("Function iterateThroughCellList wont properly iterate through CellList because of too little cells in at least one direction. \n $(Sim.CellDimensions)"; color=:red)
        
        return
    end

    ### allocate work arrays for later function calls
    Sim.CenterBox = zeros(Sim.MaxParticlesPerCell,3)
    Sim.NeighBox  = zeros(Sim.MaxParticlesPerCell,3)
    Sim.CL_Dist   = zeros(Sim.MaxParticlesPerCell^2,3)
    border=false

    for xind in x_low-1:x_high ### -1 allows do avoid double counting in subfunctions
        for yind in y_low-1:y_high
            for zind in z_low-1:z_high
                ### check whether follow up functions has to consider wrap behaviour
                border = xind==x_low || xind == x_high || yind==y_low || yind == y_high || zind==z_low || zind == z_high 
                func(Sim, xind, yind, zind,border)
                #=if border
                if xind==x_low 
                    func(Sim,x_high, yind, zind,border)
                elseif xind==x_high
                    func(Sim, x_low, yind, zind,border)
                elseif yind==y_high
                    func(Sim, xind, y_low, zind,border)                
                elseif yind==y_low
                    func(Sim,xind, y_high, zind,border)
                elseif zind==z_high
                    func(Sim, xind, yind, z_low,border)
                elseif zind==z_low
                    func(Sim,xind, yind, z_high,border)
                end
                end=#
            end
        end
    end
end

function computeDistancesForCellNeighbour(BoxLength::Vector{T}, CenterBox::Matrix{T}, NeighBox::Matrix{T}, CL_Dist::Matrix{T}, IMax::I, JMax::I, border::Bool) where {I<:Integer, T<:Real}
    ### IMax = CenterBox, JMax = NeighBox
    inv = 1.0./BoxLength

    cnt = 1 
    #=if  border
        for i in 1:IMax
            for j in 1:JMax
                for dim in 1:3
                    CL_Dist[cnt,dim] = ((NeighBox[j,dim]-CenterBox[i,dim])%BoxLength[dim])
                    #CL_Dist[cnt,dim] = CL_Dist[cnt,dim] - BoxLength[dim]*CL_Dist[cnt,dim]÷(BoxLength[dim]/2.)
                end
                CL_Dist[cnt,1] = CL_Dist[cnt,1]^2 + CL_Dist[cnt,2]^2 + CL_Dist[cnt,3]^2
                CL_Dist[cnt,1] = sqrt(CL_Dist[cnt,1])
                cnt+=1
            end
        end
    else
        for i in 1:IMax
            for j in 1:JMax
                CL_Dist[cnt,1] = (NeighBox[j,1]-CenterBox[i,1])^2
                CL_Dist[cnt,1] = (NeighBox[j,2]-CenterBox[i,2])^2 + CL_Dist[cnt,1]
                CL_Dist[cnt,1] = (NeighBox[j,3]-CenterBox[i,3])^2 + CL_Dist[cnt,1]
                CL_Dist[cnt,1] = sqrt(CL_Dist[cnt,1])
                cnt+=1
            end
        end
    end =#
    #https://en.wikipedia.org/wiki/Periodic_boundary_conditions
    for i in 1:IMax
        for j in 1:JMax
            CL_Dist[cnt,1] = abs(NeighBox[j,1]-CenterBox[i,1])
            CL_Dist[cnt,2] = abs(NeighBox[j,2]-CenterBox[i,2])
            CL_Dist[cnt,3] = abs(NeighBox[j,3]-CenterBox[i,3])

            CL_Dist[cnt,1] -= floor(I, CL_Dist[cnt,1] * inv[1]+0.5)*BoxLength[1]
            CL_Dist[cnt,2] -= floor(I, CL_Dist[cnt,2] * inv[2]+0.5)*BoxLength[2]
            CL_Dist[cnt,3] -= floor(I, CL_Dist[cnt,3] * inv[3]+0.5)*BoxLength[3]

            CL_Dist[cnt,1] =sqrt(CL_Dist[cnt,1]^2 +CL_Dist[cnt,2]^2 + CL_Dist[cnt,3]^2)
            cnt+=1
        end
    end
end