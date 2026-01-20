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


function plotContactMatrixComp(Sims::Vector{HPSAnalysis.SimData{R,I}}, Path::String, Labels::Vector{String}; Size=800,color_inter=:thermal, color_intra=:hawaii, LabelSize=14, Points=[], ranges=Vector{Vector{UnitRange{<:Integer}}}()) where {R<:Real, I<:Integer}

    comb_layout = Plots.@layout [
        a{0.975w} b{0.025w};[Plots.grid(2,2) c{0.975h, 0.025w}  ]
    ]

    InterChainMatrices = [Sim.ContactMatrices for Sim in Sims]
    addon=""
    ### replace IntraChainContactMatrix in case that 
    if length(ranges)>0
        for (i, CM) in enumerate(InterChainMatrices)
            _, InterChainMatrices[i] = compute_submatrix(CM, ranges[i])
            addon ="_with_coarse"
        end
    end

    clims_inter = (0,maximum(maximum.(InterChainMatrices)))

    intra_CMs = [sum(Sim.IntraChainContactMatrix)./length(Sim.IntraChainContactMatrix) for Sim in Sims]
    clims_intra = (minimum(minimum.(intra_CMs)), maximum(maximum.(intra_CMs)))

    data_inter = collect(range(clims_inter..., 100))
    data_intra = collect(range(clims_intra..., 100))

    fig = plot(layout=comb_layout, margins=0mm, left_margin=0mm,right_margin=0mm, size=(Size,Size))
    ### plot colorbars manually
    heatmap!(data_intra,[1],  transpose([data_intra;;]),color=color_intra,  colorbar=false, subplot=1, yaxis=nothing, xmirror=true, xlabel=L"\ln(<d_{ij}>)",xticks=true, tick_direction=:out, labelfontsize=LabelSize)
    plot!(subplot=2, framestyle=:none)
    heatmap!([1], data_inter, [data_inter;;], colorbar=false, clims=clims_inter, color=color_inter, xticks=nothing,   ymirror=true,  yguide=nothing, subplot=7, yrotation=90.0,xmirror=true, ylabel=L"P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})",xlabel=nothing, yticks=true, tick_direction=:out, right_margin=0mm, labelfontsize=LabelSize)

    Color_Mats = []
    Lims = []
    Ticks = []
    NMats = []
    for (intra_CM, inter_CM) in zip(intra_CMs, InterChainMatrices)
        matrix  = triu(inter_CM)+tril(intra_CM,-1);
        is_triu = triu(trues(size(inter_CM)))

        color_mat = [istriu ? get_color(val, clims_inter, color_inter) : get_color(val, clims_intra, color_intra)    for (istriu, val) in zip(is_triu, matrix)]
        push!(Color_Mats, color_mat)

        NMat = size(matrix,1)
        ticks = vcat([1], collect(0:50:NMat)[2:end])
        push!(NMats, NMat)
        push!(Lims, (0.5,NMat-0.5))
        push!(Ticks, ticks)
    end

    heatmap!(Color_Mats[1],yflip=false,lims=Lims[1], ticks=Ticks[1],subplot=3, xlabel=""                 , ylabel="Amino Acids j [-]" , AspectRatio=true,yrotation=90.0, tick_direction=:out)
    annotate!([NMats[1]/2.0], [NMats[1]*0.95], Labels[1], subplot=3)

    heatmap!(Color_Mats[2],yflip=false,lims=Lims[2], ticks=Ticks[2],subplot=4, xlabel=""                 , ylabel=""                  , AspectRatio=true,yrotation=90.0, tick_direction=:out)
    annotate!([NMats[2]/2.0], [NMats[2]*0.95], Labels[2], subplot=4)

    heatmap!(Color_Mats[3],yflip=false,lims=Lims[3], ticks=Ticks[3],subplot=5, xlabel="Amino Acids i [-]", ylabel="Amino Acids j [-]" , AspectRatio=true,yrotation=90.0, tick_direction=:out)
    annotate!([NMats[3]/2.0], [NMats[3]*0.95], Labels[3], subplot=5)

    heatmap!(Color_Mats[4],yflip=false, xflip=false,xlims=Lims[4],ylims=Lims[4], ticks=Ticks[4],subplot=6, xlabel="Amino Acids i [-]", ylabel=""                  , AspectRatio=true,yrotation=90.0,  tick_direction=:out)
    annotate!([NMats[4]/2.0], [NMats[4]*0.95], Labels[4], subplot=6, legend=nothing)

    cnt = 8
    for (i,j) in enumerate(3:6)
        if length(Points)>0 && length(Points[i])>0
            plot!( twinx(fig[j]);) # creates a right‑hand axis that shares the same limits
            plot!( twiny(fig[j]);) # creates a upper‑hand axis that shares the same limits

            plot!(yticks=(Points[i], ["" for _ in Points[4]]), yforeground_color_axis=:red, subplot=cnt  , lims=Lims[i], tick_direction=:out,  foreground_color_border=:white) 
            plot!(xticks=(Points[i], ["" for _ in Points[4]]), xforeground_color_axis=:red, subplot=cnt+1, lims=Lims[i],tick_direction=:out,  foreground_color_border=:white) 
            cnt += 2
        end
    end

    Plots.savefig(fig, Path*"Comp_ContactMatrix$(addon).png")
    Plots.savefig(fig, Path*"Comp_ContactMatrix$(addon).pdf")
    return fig
end
