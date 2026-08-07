using LaTeXStrings, Plots

function plotBondAngleHistComp(Sims::Vector{HPSAnalysis.SimData{R,I}}, Path::String) where {R <: Real , I<: Integer}
    BondAnglePlot = Plots.plot(xlabel="Bond Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Bond angle histogram")
    for Sim in Sims
        Plots.plot!(axes(Sim.BondAngleHist,1)*180.0/size(Sim.BondAngleHist,1), Sim.BondAngleHist[:]*100, label=Sim.SimulationName)
    end
    Plots.savefig(BondAnglePlot, Path*"Comp_BondAngleHist.png")
    Plots.savefig(BondAnglePlot, Path*"Comp_BondAngleHist.pdf")
    return BondAnglePlot
end

function plotRGAutocorrComp(Sims::Vector{HPSAnalysis.SimData{R,I}},  Path::String) where {R<:Real, I<:Integer}
    RG_Auto= Plots.plot(xlabel="lag time [τ]", ylabel="autocorr []", xlims=(0,500),ylims=(-0.5, 1.0))
    for Sim in Sims
        axis = axes(Sim.RGAutocorr[1,:])
        avg = sum(Sim.RGAutocorr, dims=1)[1,:]./Sim.NChains
        err = zeros_like(avg)
        for chain in 1:Sim.NChains
            err .+= (Sim.RGAutocorr[chain,:] .- avg).^2
        end
        err .= sqrt.(err./Sim.NChains)./(Sim.NChains-1)

        Plots.plot!(axis, avg,yerr=err,  label=Sim.SimulationName)
    end

    Plots.savefig(RG_Auto, Path*"_RGAutocorrComp.png")
    Plots.savefig(RG_Auto, Path*"_CompRGAutocorrComp.pdf")
    return RG_Auto
end

function plotDihedralAngleHistComp(Sims::Vector{HPSAnalysis.SimData{R,I}}, Path::String) where {R<:Real, I<:Integer}
    DihedralAnglePlot= Plots.plot(xlabel="Diehedral Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Dihedral angle histogram", label="",legend_position=:topleft)
    for Sim in Sims
        DihedralAnglePlot= Plots.plot!(axes(Sim.TorsionHist,2), Sim.TorsionHist[end,:]*100,  label=Sim.SimulationName)
    end
    #Plots.plot!(( x)*180.0/π, (P+P2)*100/2.0, label="Probability")
    Plots.savefig(DihedralAnglePlot, Path*"Comp_DihedralAngleHist.png")
    Plots.savefig(DihedralAnglePlot, Path*"Comp_DihedralAngleHist.pdf")
    return DihedralAnglePlot
end

function computeContactMatrixColorMatrix(Sims::Vector{HPSAnalysis.SimData{R,I}},ranges::Vector{Vector{UnitRange{I2}}}=Vector{Vector{UnitRange{Int32}}}(),color_inter=:thermal,color_intra=:hawaii, xspace::Int64=50) where {R<:Real, I<:Integer,I2<:Integer}
   # Prepare the contact matrices (unchanged from the original code)
    inter_CMs = [Sim.ContactMatrices for Sim in Sims]
    addon = ""

    if length(ranges) > 0
        for (i, CM) in enumerate(inter_CMs)
            _, inter_CMs[i] = compute_submatrix(CM, ranges[i])
            addon = "_with_coarse"
        end
    end

    clims_inter = (0, maximum(maximum.(inter_CMs)))
    intra_CMs   = [sum(Sim.IntraChainContactMatrix) ./ length(Sim.IntraChainContactMatrix) for Sim in Sims]
    clims_intra = (minimum(minimum.(intra_CMs)), maximum(maximum.(intra_CMs)))

    data_inter = collect(range(clims_inter..., 100))
    data_intra = collect(range(clims_intra..., 100))

    # Build colour matrices, limits, ticks for *each* simulation
    Color_Mats = Vector{Matrix{RGB{Float64}}}()
    Lims       = Vector{Tuple{Float64,Float64}}()
    Ticks      = Vector{Vector{Int}}()
    NMats      = Vector{Int}()

    for (intra_CM, inter_CM) in zip(intra_CMs, inter_CMs)
        matrix = triu(inter_CM) + tril(intra_CM, -1)
        is_triu = triu(trues(size(inter_CM)))
        colour_mat = [istriu_ ? get_color(val, clims_inter, color_inter) : get_color(val, clims_intra, color_intra) for (istriu_, val) in zip(is_triu, matrix) ]

        push!(Color_Mats, colour_mat)
        NMat = size(matrix, 1)
        push!(NMats, NMat)
        push!(Lims, (0.5, NMat - 0.5))
        push!(Ticks, vcat([1], collect(0:xspace:NMat)[2:end]))
    end

    return inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks, addon

end

function plotContactMatrixComp(Sims::Vector{HPSAnalysis.SimData{R,I}},Path::String,Labels::Vector{String};Size=(600,600),color_inter=:thermal,color_intra=:hawaii,LabelSize= 10,Points::Vector{Vector{R2}}=Vector{Vector{Int32}}(),xspace::Int64=50,ranges::Vector{Vector{UnitRange{I2}}}=Vector{Vector{UnitRange{Int32}}}(),Grid=nothing,topsepspace=-7mm, bspace=0mm,rightsepspace=-7mm, tspace=0mm,lspace=0mm,inPercent=false) where {R<:Real, I<:Integer,I2<:Integer, R2<:Real}

    nSims = length(Sims) # how many simulations we have

    # Determine the grid that will hold the heat‑maps
    ncols = isnothing(Grid) ? ceil(Int, sqrt(nSims)) : Grid[1]
    nrows = isnothing(Grid) ? ceil(Int, nSims/ncols) : Grid[2] 
    grid = Plots.grid(nrows,ncols)
    comb_layout = @layout [
        a{0.975w}  b{0.025w};
        [grid]    c{0.975h,0.025w}
    ]
    nGrid=ncols*nrows

    inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks,addon = computeContactMatrixColorMatrix(Sims,ranges,color_inter,color_intra, xspace)

    # Figure and colour‑bars
    fig = plot(layout = comb_layout,margins = -3mm, left_margin = 0mm, right_margin = 0mm,size = Size, top_margin=0mm, bottom_margin=0mm, minorticks=10)

    # intra‑chain colour bar (left)
    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false,subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\ln(\langle d_{ij}\rangle/\AA)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=topsepspace, clims=clims_intra,top_margin=tspace)

    # empty placeholder (keeps the left bar thin)
    plot!(subplot = 2, framestyle = :none, bottom_margin=topsepspace, left_margin=rightsepspace,right_margin=2mm)

    ylabel_ = inPercent ?  L"P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})\quad [\%]" :  L"P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})" 
    # inter‑chain colour bar (right)
    clims_inter = inPercent ? clims_inter.*100.0 : clims_inter
    data_inter .*= inPercent ? 100.0 : 1.0
    heatmap!([1], data_inter, [data_inter;;],colorbar = false, clims = clims_inter,color = color_inter,xticks = nothing,ymirror = true,yguide = nothing,subplot = nGrid + 3,yrotation = 90.0,xmirror = true,ylabel = ylabel_,xlabel = nothing,yticks = true,tick_direction = :out,right_margin = 2.5mm,labelfontsize = LabelSize, left_margin=rightsepspace,yminorticks=10)

    # Plot the heat‑maps - one subplot per simulation
    for i in 1:nSims
        subplot_idx = 2 + i                      # 3,4,5,… (first heat‑map is 3)

        # Decide whether we need axis‑labels for this cell
        row = ((i - 1) ÷ ncols) + 1
        col = ((i - 1) % ncols) + 1

        xlbl = (row == nrows) ? "Amino Acids i " : ""
        ylbl = (col == 1)     ? "Amino Acids j " : ""
        bs = (row == nrows) ? bspace  : 0mm
        ls = (col ==1) ? lspace : 0mm

        heatmap!(Color_Mats[i],yflip = false,lims = Lims[i],ticks = Ticks[i],subplot = subplot_idx, xlabel = xlbl, ylabel = ylbl, AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top,bottom_margin=bs,left_margin=ls)

        # label in the upper‑right corner of each panel
        annotate!([NMats[i] / 2.0],[NMats[i] * 0.95],text(get(Labels, i, ""), :black, :center, LabelSize), subplot = subplot_idx)
    end

    for i in 3+nSims:2+nGrid
        row = ((i - 1) ÷ ncols) + 1
        col = ((i - 1) % ncols) + 1
        xlbl = (row == nrows) ? "Amino Acids i " : ""
        ylbl = (col == 1)     ? "Amino Acids j " : ""
        plot!(subplot=i, xlabel = xlbl,framestyle = :none,)
    end

    # (Optional) plot point‑sets on top of each heat‑map
    if length(Points) > 0
        cnt = nSims + 4
        for (i, j) in enumerate(3:(nSims + 2))   # `j` is the heat‑map subplot index
            len = NMats[i]
            if length(Points) > 0 && length(Points[i]) > 0
                # Create twin axes that share the limits of the heat‑map
                for point in Points[i]
                    plot!([point,point],[floor(I,len*0.975), len], label="", c=:red, subplot=j, linewidth=2.0)
                    plot!([floor(I,len*0.975), len],[point,point], label="", c=:red, subplot=j, linewidth=2.0)
                end
            end
        end
    end

    Plots.savefig(fig, Path * "Comp_ContactMatrix$(addon).png")
    Plots.savefig(fig, Path * "Comp_ContactMatrix$(addon).pdf")
    return fig
end


function plotContactMatrixDifference(Sims::Vector{HPSAnalysis.SimData{R,I}}, filename::String;color_inter=:thermal,color_intra=:hawaii,ranges::Vector{UnitRange{I2}}=Vector{UnitRange{I2}}(),xspace::Int64=50,Size=(1050,400), LabelSize=10, bspace=-3mm,lspace=-6mm, rspace=-3mm, tspace=7mm, Labels::Vector{String}=["",""],CheckSequence=true) where {R<:Real, I<:Integer,I2<:Integer}
    @assert length(Sims)>=2 "Number of Simulations needs to be atleast 2"

     inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks,addon = computeContactMatrixColorMatrix(Sims,Vector{Vector{UnitRange{Int32}}}(),color_inter,color_intra, xspace)
    grid = Plots.grid(1,3)

    comb_layout = @layout [
       [a{0.5w} b{0.5w,0.06h};  
       e{0.5w} f{0.5w,0.9h}];
    ]

    ### SUBPLOT 1

    fig1 = Plots.plot(size=Size,bottom_margin=5mm, layout=comb_layout)

    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false,subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\ln(<d_{ij}>)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=bspace, clims=clims_intra, top_margin=tspace)


    heatmap!(data_inter, [1],transpose([data_inter;;]),color = color_inter,colorbar = false,subplot = 2,yaxis = nothing,xmirror = true,xlabel = L"P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=bspace, clims = clims_inter, top_margin=tspace)

    xlbl = "Amino Acids i "
    ylbl = "Amino Acids j "

    Plots.heatmap!(Color_Mats[1], subplot=3,yflip = false,lims = Lims[1],ticks = Ticks[1], xlabel = xlbl, ylabel = ylbl, AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top,left_margin=4mm)

    Plots.heatmap!(Color_Mats[2], subplot=4, yflip = false,lims = Lims[2],ticks = Ticks[2], xlabel = xlbl, ylabel = "", AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top)

    # label in the upper‑right corner of each panel
    annotate!([NMats[1] / 2.0],[NMats[1] * 0.95],text(get(Labels, 1, ""), :black, :center, LabelSize), subplot = 3)
    annotate!([NMats[2] / 2.0],[NMats[2] * 0.95],text(get(Labels, 2, ""), :black, :center, LabelSize), subplot = 4)

    ### SUBPLOT 2

    comb_layout = @layout [[ c{0.95w} d{0.06h,0.025w};
                             g{0.95w} h{0.9h ,0.025w}];]
    fig2=Plots.plot(size=Size,bottom_margin=5mm, layout=comb_layout)#

    ranges = length(ranges) > 0 ? ranges :  [1:size(inter_CMs[1])[1], 1:size(inter_CMs[2])[1]]
    @assert length(ranges[1]) == length(ranges[2]) "Matrix areas specified with arguemtns \"ranges\" need to be of same size"
    if CheckSequence
        @assert Sims[1].Sequences[1][ranges[1]]==Sims[2].Sequences[1][ranges[2]] "Ranges do not contain the same sequence"
    end

    inter_A = inter_CMs[1][ranges[1],ranges[1]]
    inter_B = inter_CMs[2][ranges[2],ranges[2]]

    intra_A = intra_CMs[1][ranges[1],ranges[1]]
    intra_B = intra_CMs[2][ranges[2],ranges[2]]

    intra_diff = intra_A.- intra_B 
    inter_diff = inter_A.- inter_B 

    matrix = triu(inter_diff) + tril(intra_diff, -1)
    is_triu = triu(trues(size(inter_diff)))

    cmax =maximum(maximum(abs.(inter_diff)))
    clims_inter = (-cmax, cmax)
    cmax =maximum(maximum(abs.(intra_diff)))
    clims_intra = (-cmax, cmax)

    color_inter = :romaO
    color_intra = :tarn

    colour_mat = [istriu_ ? get_color(val, clims_inter, color_inter) : get_color(val, clims_intra, color_intra) for (istriu_, val) in zip(is_triu, matrix) ]
    
    ticks = vcat([1], collect(0:xspace:size(inter_diff)[1])[2:end])
    lims= (1,size(inter_diff)[1])

    data_inter = collect(range(clims_inter..., 100))
    data_intra = collect(range(clims_intra..., 100))

    rr_space=-5mm
    sep_space=-3mm
    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false, subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\Delta\ln(<d_{ij}>)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=bspace, clims=clims_intra, top_margin=tspace,left_margin=sep_space,right_margin=rr_space)

    Plots.heatmap!(colour_mat, subplot=3, yflip = false, ticks=ticks, xlabel = xlbl, ylabel = "", AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top,lims=lims, left_margin=sep_space,right_margin=rr_space)



    # empty placeholder (keeps the right bar thin)
    plot!(subplot = 2, framestyle = :none, bottom_margin=bspace, left_margin=lspace,right_margin=rspace, top_margin=tspace)


    heatmap!([1], data_inter, [data_inter;;],colorbar = false, clims = clims_inter,color = color_inter,xticks = nothing,ymirror = true,yguide = nothing,subplot = 4,yrotation = 90.0,xmirror = true,ylabel = L"\Delta P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})",xlabel = nothing,tick_direction = :out,labelfontsize = LabelSize, left_margin=lspace,yminorticks=10,right_margin=rspace,)

    for (i,range) in zip([3,4],ranges)
        low  = first(range)
        high =  last(range)
        plot!(fig1,[low, low,high, high,low], [low, high, high, low,low], linewidth=3, linecolor=:red, label="", linestyle=:solid, subplot=i)
    end

    ### COMBINE PLOTS
    fig3 =vline([1], linewidth=10,label="", c=:black, framestyle=:none,left_margin=-9mm)
    comb_layout= @layout [ j{0.63w}  b{0.025w} c{0.3w,1.0h}]

    fig = Plots.plot(fig1,fig3,fig2,layout=comb_layout)

    Plots.savefig(fig, filename*".png")
    Plots.savefig(fig, filename*".pdf")

    return fig
end
