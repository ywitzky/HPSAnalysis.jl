
#module GSD
#include("./Unify.jl")

include("./Plot/Comp_Plots.jl")

using Plots, Printf, Measures
#

#=
function plotRGAutocorr(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    x = axes(Sim.RGAutocorr,2)
    avg= reduce(+,Sim.RGAutocorr, dims=1)[1,:]./Sim.NChains

    Plots.plot(xlabel="Lag time [steps]", ylabel="AutoCorr [steps]")
    for chain in 1:Sim.NChains
        Plots.plot!(x, Sim.RGAutocorr[chain,:], label="")#, label="$(chain)")
    end
    Plots.vline!([Sim.RGMeasureStep], label="")
    Plots.hline!([0.0], label="")


    Plots.plot!(x, avg, color=:black, label="avg")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.pdf")
end
=#

function plotREEVecAutocorr(Sim::LammpsAnalysis.SimData{R,I}, conf=0.1) where {R<:Real, I<:Integer}
    x = (collect(R, axes(Sim.REEVecSeries,1)).-convert(R,1.0)) * convert(R,conf)

    Plots.plot(xlabel="Lag time [ns]", ylabel="AutoCorr [steps]")
    avg = sum(Sim.REEVecSeries[:, :, 1], dims=2)./Sim.NChains
    for chain in 1:Sim.NChains
        Plots.plot!(x, Sim.REEVecSeries[:, chain, 1], label="")
    end
    Plots.plot!(x,avg[:, 1, 1], label="avg.", c=:black)

    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr.pdf")
    Plots.plot!(xlim=(0,3.0))
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr_zoom.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr_zoom.pdf")
end

#=function plotREEAutocorr(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    x = axes(Sim.REEAutocorr,2)
    avg= reduce(+,Sim.REEAutocorr, dims=1)[1,:]./Sim.NChains

    Plots.plot(xlabel="Lag time [steps]", ylabel="AutoCorr [steps]")
    for chain in 1:Sim.NChains
        Plots.plot!(x, Sim.REEAutocorr[chain,:], label="")#, label="$(chain)")
    end
    Plots.vline!([Sim.RGMeasureStep], label="")
    Plots.hline!([0.0], label="")

    Plots.plot!(x, avg, color=:black, label="avg")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.pdf")
end=#

function plotChargeBondTime(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.CellResolution=20.
    Sim.ChargeAnalysisCutoff=20.
    if length(Sim.ChargeContactTimeHist)==0
        SetFieldsFromData(Sim, ["ChargeContactTimeHist", "NegChargeToHistID", "PosChargeToHistID"])
       # computeChargeCorrelations(Sim)
    end

    #
   # Plots.plotlyjs()

    ChargeBondTime =  Plots.plot(size = (1150, 500),  yticks=([10^-8,10^-7,10^-6,10^-5,10^-4,0.001,0.01,0.1,1,10,100], ["  10^-8","  10^-7","  10^-6","  10^-5","  10^-4","  10^-3","  10^-2","  10^-1","  10^0","  10^1","  10^2"]),yaxis=:log, xaxis=:log, xlabel="Time ns", ylabel="P(r<2nm) [%]", grid=true, gridalpha=0.25, legend=:topright)

	step_to_time_in_ns = 10*100_000*10^-15/10^-9
	range = (axes(Sim.ChargeContactTimeHist,3)*step_to_time_in_ns)
	time = collect(0:1.0:length(axes(Sim.ChargeContactTimeHist,3))-1)

    for pos in keys(Sim.PosChargeToHistID)
		pos_id = Sim.PosChargeToHistID[pos]
		for neg in keys(Sim.NegChargeToHistID)
			neg_id = Sim.NegChargeToHistID[neg]
			hist = Sim.ChargeContactTimeHist[pos_id, neg_id,:]/sum(Sim.ChargeContactTimeHist[pos_id, neg_id,:])*100
			time_avg = sum(time.*hist)*step_to_time_in_ns/100
            Plots.plot!(range[1:150],hist[1:150], label="Charges: $(pos) , $(neg)  \n avg time: $(@sprintf("%.3f", time_avg)) ns", markershape=:cross )
		end
	end
    Plots.savefig(ChargeBondTime, Sim.PlotPath*Sim.SimulationName*"_ChargeBonds.png")

    Plots.savefig(ChargeBondTime, Sim.PlotPath*Sim.SimulationName*"_ChargeBonds.pdf")

	return ChargeBondTime
end

function plotChargeBonds(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    if length(Sim.ChargeChainContactMatrix)==0
        SetFieldsFromData(Sim, ["ChargeChainContactMatrix", "MaxChainLength", "ChargeResidueContactMatrix"])
    #    computeChargeCorrelations(Sim)
    end
    

    Plots.gr() ### change Plots backend to get indepenent colorbars
	ChainMatrix = Plots.heatmap(1:Sim.NChains,1:Sim.NChains,Sim.ChargeChainContactMatrix, xlabel="Chain i", ylabel="Chain j" , title="Chain Contacts", colorbar=true)
	ResidueMatrix = Plots.heatmap(1:Sim.MaxChainLength,1:Sim.MaxChainLength,Sim.ChargeResidueContactMatrix, xlabel="Residue i", ylabel="Residue j" , title="Residue Contacts", colorbar=true)
	
	ChargeMatrix  = Plots.plot(ChainMatrix, ResidueMatrix, layout=2, size=(1200,500), left_margin=8mm, bottom_margin=5mm)

    Plots.savefig(ChargeMatrix, Sim.PlotPath*Sim.SimulationName*"_ChargeBondMatrix.png")
    Plots.savefig(ChargeMatrix, Sim.PlotPath*Sim.SimulationName*"_ChargeBondMatrix.pdf")

    return ChargeMatrix
end

function plotMSDofChains(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Plots.gr()
    MSD_Plot=Plots.plot()
    MSD_Diff=Plots.plot()


    #MinLength = minimum(Sim.BoxLength)/2. ### visualisation of a improper behaviour of the Analysis after certain diffusion values
    #MaxDeviation = MinLength^2
    MSD_all = Plots.plot()
    for chain in 1:Sim.NChains			
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,1]+data[:,2]+data[:,3], label="") 
    end


    N = length(axes(Sim.MSD,1))
    avg = sum(Sim.MSD, dims=1)[1,:,:]/Sim.NChains
    steps = length(axes(avg,1))
    diff_avg  = avg[2:steps,:]-avg[1:steps-1,:]

    avg_err =zero(avg) ### works like zeros_like in python
    diff_avg_err =zero(diff_avg) ### works like zeros_like in python
    diff = Sim.MSD[:,2:steps,:]-Sim.MSD[:,1:steps-1,:]
    for j in axes(Sim.MSD,3)
        for i in axes(Sim.MSD,2)
            avg_err[i,j] =  sum((avg[i,j] .- Sim.MSD[:,i,j]).^2)
        end
        for i in axes(diff,2)
            diff_avg_err[i,j] =  sum((diff_avg[i,j] .- diff[:,i,j]).^2)
        end
    end
    avg_err .= sqrt.(avg_err./N)/sqrt(N-1)
    diff_avg_err .= sqrt.(diff_avg_err./N)/sqrt(N-1)

    MaxVal = maximum(avg[1: ceil(I,Sim.NSteps*0.75)])
    ylims = [0, MaxVal]

    width=avg_err[:,1].+avg_err[:,2].+avg_err[:,3]
    Plots.plot!(MSD_all,axes(avg,1), avg[:,1]+avg[:,2]+avg[:,3]+width,fillrange=avg[:,1]+avg[:,2]+avg[:,3]-width, color=:gray, label="uncertainty")
    Plots.plot!(MSD_all,axes(avg,1), avg[:,1]+avg[:,2]+avg[:,3], color=:black, label="avg.", title="X&Y&Z-component", ylim=ylims.*3, legend=:bottomright)
    

    Plots.plot()
    for chain in 1:Sim.NChains
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,1], label="") #+data[:,2]
    end
    Plots.plot!(axes(avg,1), avg[:,1]+avg_err[:,1],fillrange=avg[:,1]-avg_err[:,1], color=:gray, label="")
    MSD_XY = Plots.plot!( axes(avg,1), avg[:,1], color=:black, title="X-component" , label="", ylim=ylims) #+avg[:,2], title="X-Y"
    

    Plots.plot()
    for chain in 1:Sim.NChains
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,2], label="") #+data[:,3]
    end
    Plots.plot!(axes(avg,1), avg[:,2]+avg_err[:,2],fillrange=avg[:,2]-avg_err[:,2], color=:gray, label="")
    MSD_XZ = Plots.plot!( axes(avg,1), avg[:,2], color=:black, title="Y-component" , label="", ylim=ylims)

    
    Plots.plot()
    for chain in 1:Sim.NChains
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,3], label="") #data[:,2]+
    end
    Plots.plot!(axes(avg,1), avg[:,3]+avg_err[:,3],fillrange=avg[:,3]-avg_err[:,3], color=:gray, label="")
    MSD_YZ = Plots.plot!(axes(avg,1), avg[:,3], color=:black, title="X-component" , label="" , ylim=ylims)
    MSD_Plot = Plots.plot(MSD_all, MSD_XY,MSD_XZ,MSD_YZ, layout=4, plot_title="MSD") # xscale=:log, yscale=:log,


 
    #max_step=Int32(min(500., Float32(length(axes(avg,1)))))
    #lims = [minimum(diff_avg[1:max_step,:])*0.95, maximum(diff_avg[1:max_step,:])*1.05]
    diff_ylims = [0, maximum(diff_avg[1: ceil(I,Sim.NSteps*0.75)])]
    Diff_all = Plots.plot()
    
    val=diff_avg[:,1]+diff_avg[:,2]+diff_avg[:,3]
    width =diff_avg_err[:,1]+diff_avg_err[:,2]+diff_avg_err[:,3]
    Diff_all = Plots.plot!(1:steps-1, val+width, fillrange=val-width, color=:gray, label="")
    Plots.plot!(1:steps-1, val, color=:black, label="" , title="X&Y&Z-component", ylims= diff_ylims.*3)


    Diff_XY  = Plots.plot(1:steps-1, diff_avg[:,1]+diff_avg_err[:,1], fillrange=diff_avg[:,1]-diff_avg_err[:,1], color=:gray, label="" , title="X-component", ylims= diff_ylims)
    Plots.plot!(1:steps-1, diff_avg[:,1], color=:black, label="" , title="X", ylims= diff_ylims)

    Diff_XZ  = Plots.plot(1:steps-1, diff_avg[:,2]+diff_avg_err[:,2], fillrange=diff_avg[:,2]-diff_avg_err[:,2], color=:gray, label="" , title="Y-component", ylims= diff_ylims)
    Plots.plot!(1:steps-1, diff_avg[:,2], color=:black, label="" , title="Y", ylims= diff_ylims)

    Diff_YZ  = Plots.plot(1:steps-1, diff_avg[:,3]+diff_avg_err[:,3], fillrange=diff_avg[:,3]-diff_avg_err[:,3], color=:gray, label="" , title="Z-component", ylims= diff_ylims)
    Plots.plot!(1:steps-1, diff_avg[:,3], color=:black, label="" , title="Z", ylims= diff_ylims)
    annotate!(Diff_all,(-0.2, 1.3), text("$(Sim.SimulationName)")) ### hacky way to get plot title that isnt set for subplots

    
    MSD_Diff = Plots.plot(Diff_all, Diff_XY, Diff_XZ, Diff_YZ, layout=4, plot_title="d/dt MSD")
    MSD_Master = Plots.plot(MSD_Plot , MSD_Diff, layout=2, size=(1600, 600))#,yscale=:log,xscale=:log)

    Plots.savefig(MSD_Master, Sim.PlotPath*Sim.SimulationName*"_MSD_Chains.png")
    Plots.savefig(MSD_Master, Sim.PlotPath*Sim.SimulationName*"_MSD_Chains.pdf")
    return MSD_Master
end

function plotBondHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

    #max_val = maximum(axes(Sim.BondHist),1)
    BondPlot= Plots.plot(axes(Sim.BondHist,1)/100, Sim.BondHist*100, xlabel="Bond distance r [AA]", ylabel="P(r) [%]", label="", xlim=(3,4.5))

    Plots.savefig(BondPlot, Sim.PlotPath*Sim.SimulationName*"_BondHist.png")
    Plots.savefig(BondPlot, Sim.PlotPath*Sim.SimulationName*"_BondHist.pdf")
    return BondPlot
end

function plotDihedralAngleHist(Sim)
    kα1 = 11.4
    kα2 = 0.15
    θα1 = 0.90
    θα2 = 1.02
    kβ1 = 1.80
    kβ2 = 0.65
    θβ1 = -1.55
    θβ2 = -2.50
    e0  =0.27
    e1  = 0.14
    e2  = 0.40
    Dihed_a(x, ed) = @. exp(-kα1*(x-θα1)^2-ed) + exp(-kα2*(x-θα2)^4+e0) + exp(-kα2*(x-θα2+2*π)^4+e0)
    Dihed_b(x, ed) = @. exp(-kβ1*(x-θβ1)^2+e1+ed) + exp(-kβ1*(x-θβ1-2*π)^2+e1+ed) + exp(-kβ2*(x-θβ2)^4+e2) + exp(-kβ2*(x-θβ2-2*π)^4+e2)
    U(x,ed) = @. -log(Dihed_a(x, ed)+Dihed_b(x, ed))

    NA = 6.022
    kb = 1.38
    kcal = 4184
    kt = kb *Sim.TargetTemp *NA/kcal

    x = LinRange(-π,π,500)
    E = U(x, -1.0)
    W = exp.(-1.0*E/kt)#.*sin.( x)
    P = W./sum(W)*(2*π/500.0)^-1

    E = U(x, 1.0)
    W = exp.(-1.0*E/kt)#.*sin.( x)
    P2 = W./sum(W)*(2*π/500.0)^-1


    DihedralAnglePlot= Plots.plot(axes(Sim.TorsionHist,2), Sim.TorsionHist[end,:]*100, xlabel="Diehedral Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Dihedral angle histogram", label="")
    Plots.plot!(( x)*180.0/π, (P+P2)*100/2.0, label="Probability")

    #DihedralAnglePlot= Plots.plot( xlabel="Diehedral Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Dihedral angle histogram", label="")
    #=for i in 2:length(axes(Sim.TorsionHist,1))
        Plots.plot!(axes(Sim.TorsionHist,2), Sim.TorsionHist[i,:]*100, label="")
    end=#
    Plots.savefig(DihedralAnglePlot, Sim.PlotPath*Sim.SimulationName*"_DihedralAngleHist.png")
    Plots.savefig(DihedralAnglePlot, Sim.PlotPath*Sim.SimulationName*"_DihedralAngleHist.pdf")
    return DihedralAnglePlot
end

function plotBondAngleHist(Sim)
    BondAnglePlot= Plots.plot(axes(Sim.BondAngleHist,1)*180.0/size(Sim.BondAngleHist,1), Sim.BondAngleHist[:]*100, xlabel="Bond Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Bond angle histogram", label="")

    Plots.savefig(BondAnglePlot, Sim.PlotPath*Sim.SimulationName*"_BondAngleHist.png")
    Plots.savefig(BondAnglePlot, Sim.PlotPath*Sim.SimulationName*"_BondAngleHist.pdf")
    return BondAnglePlot
end

function plotRGHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Plots.plot()
    RG_Hist = Plots.histogram([(Sim.RGSeries...)...], normalize=:pdf, xlabel="RG [AA]", ylabel="P(RG)", label="") #,, bins=b_range color=:gray)
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGHist.png")
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGHist.pdf")
    return RG_Hist
end

function plotRGSeries(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    RG_Hist= Plots.plot()
    axis = axes(Sim.RGSeries[1,:])
    for chain in 1:Sim.NChains
        RG_Hist = Plots.plot!(axis, Sim.RGSeries[chain,:], xlabel="t [?]", ylabel="RG [AA]", label="")
    end
    Plots.vline!([Sim.EquilibrationTime], label="Choosen EquilibrationTime")
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGSeries.png")
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGSeries.pdf")
    return RG_Hist
end

function plotRGAutocorr(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    RG_Auto= Plots.plot(xlabel="lag time [τ]", ylabel="autocorr []", xlims=(0,500),ylims=(-0.5, 1.0))
    axis = axes(Sim.RGAutocorr[1,:])
    avg = sum(Sim.RGAutocorr, dims=1)[1,:]./Sim.NChains

    for chain in 1:Sim.NChains
        RG_Auto = Plots.plot!(axis, Sim.RGAutocorr[chain,:], label="")
    end
    Plots.plot!(axis, avg,  label="avg. over chains", c=:black)

    #Plots.vline!([Sim.EquilibrationTime], label="Choosen EquilibrationTime")
    Plots.savefig(RG_Auto, Sim.PlotPath*Sim.SimulationName*"_RGAutocorr.png")
    Plots.savefig(RG_Auto, Sim.PlotPath*Sim.SimulationName*"_RGAutocorr.pdf")
    return RG_Auto
end

function plotREEHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Plots.plot()
    REE_Hist = Plots.histogram([(Sim.REESeries...)...], normalize=:pdf, xlabel="REE [AA]", ylabel="P(REE)", label="") #,, bins=b_range color=:gray)
    Plots.savefig(REE_Hist, Sim.PlotPath*Sim.SimulationName*"_REEHist.png")
    Plots.savefig(REE_Hist, Sim.PlotPath*Sim.SimulationName*"_REEHist.pdf")
    return REE_Hist
end

function animateDensityGif(Sim::SimData{R,I}, NBins=100) where {R<:Real, I<:Integer}

    conversion = 1.66053906660/(Sim.BoxLength[1]*Sim.BoxLength[2]/Float32(NBins)^2*Sim.BoxLength[3])
    weights= Sim.Masses*conversion
    anim = @animate for step ∈ 1:10:Sim.NSteps
        fig=Plots.plot(layout=@layout[a{0.35w} b{0.2w} c{0.45w}] , size=(1750/1.6,500/2), left_margin=5mm, right_margin=2.5mm, margin=0mm, bottom_margin=6mm)
        Plots.histogram2d!(Sim.y[:,step], Sim.x[:, step], weights=weights, nbins = NBins, show_empty_bins=true, clims=(0,1.5),  xlabel="y", ylabel="x", subplot=1, colorbar=:false,xlims=Tuple(Sim.BoxSize[2,:]), ylims=Tuple(Sim.BoxSize[1,:]))
        Plots.histogram2d!(Sim.x[:, step], Sim.z[:,step], weights=weights, nbins = NBins, show_empty_bins=true, clims=(0,1.5), xlabel="x", ylabel="z", subplot=2, colorbar=:false, xlims=Tuple(Sim.BoxSize[1,:]), ylims=Tuple(Sim.BoxSize[3,:]))
        Plots.histogram2d!(Sim.y[:, step], Sim.z[:,step], weights=weights, nbins = NBins, show_empty_bins=true, clims=(0,1.5), xlabel="y", ylabel="z", subplot=3, colorbar=:right,colorbar_title="kg/L", xlims=Tuple(Sim.BoxSize[2,:]), ylims=Tuple(Sim.BoxSize[3,:]))
        Plots.annotate!(1.55*Sim.BoxSize[2,2], 1.15*Sim.BoxSize[3,1], "$(step/1000.0) μs", subplot=3)
        fig
    end
    gif(anim, "/localscratch/HPS_DATA/HPS-Alpha/SLAB/DensityGifs/$(Sim.SimulationName)_$(Sim.TargetTemp).gif", fps = 30)
end

function plotAvgSlabDensity(Sim::SimData{R,I}; Windowlength=100) where {R<:Real, I<:Integer}
    if Windowlength > Sim.NSteps
        Windowlength=Sim.NSteps
    end

    #local xticks=([-150,-100,-50,0,50,100,150], ["-750", "-500", "-250", "0", "250", "500","750"])

    xaxis = axes(Sim.SlabHistogramSeries)[1]
    ### newer iterations dont measure at every frame. Detect which frames actually contain data
    NMeasurements= sum(Sim.SlabHistogramSeries[0,Sim.NSteps-Windowlength:Sim.NSteps,1].!=0.0)
    println("NMeas: $NMeasurements")
    AvgHist = sum(Sim.SlabHistogramSeries[xaxis,Sim.NSteps-Windowlength:Sim.NSteps,:], dims=2)./(NMeasurements)

    #layout=(length(Temperatures),length(Proteins)  )
    # xticks=xticks,, ylim=(0,0.55)
    fig = Plots.plot(dpi=300, ylabel= "avg. density"* "  [kg/L]" , xlabel= "z-Axis [Å]" ) #,size=(600,300), figsize=(12cm, 4.5cm)
    Plots.plot!(xaxis, AvgHist[:,1,1], color=:black, label="",  title=Sim.SimulationName)# ind ==1 ? Prot : "")#, bottom_margin= ind==3 ? 10mm : 0mm, left_margin= Num==1 ? 15mm : 0mm ) #title= Num ==1 ? Prot : "phosphorylated RS41")

    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_AvgSlabHist_LastFrames.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_AvgSlabHist_LastFrames.pdf")
    return fig
end

function plotAvgSlabDensityEvolution(Sim::SimData{R,I}; Windowlength=100) where {R<:Real, I<:Integer}
    if Windowlength > Sim.NSteps
        Windowlength=Sim.NSteps
    end
    
    xval= LinRange(Sim.BoxSize[Sim.SlabAxis,1], Sim.BoxSize[Sim.SlabAxis,2], 9)
    xticks = (xval.-Sim.BoxSize[Sim.SlabAxis,1], string.(xval))
    xaxis = axes(Sim.SlabHistogramSeries)[1]

    tmp = ceil(I,Sim.NSteps/Windowlength)
    AvgHist = zeros(R,(size(Sim.SlabHistogramSeries)[1], size(Sim.SlabHistogramSeries)[3], tmp) )
    for i in 1:tmp-1
        start = 1 + (i-1)*Windowlength
        stop = start + Windowlength-1
        ### newer iterations dont measure at every frame. Detect which frames actually contain data
        NMeasurements= sum(Sim.SlabHistogramSeries[0,start:stop,1].!=0.0)
        AvgHist[:,:,i] .= OffsetArrays.no_offset_view((reduce(+, Sim.SlabHistogramSeries[xaxis,start:stop,:], dims=2)./(NMeasurements)))[:,1,:] ### convert offsetarray back so that we can assign it to AvgHist
    end
    start = (tmp-1)*Windowlength+1
    stop=Sim.NSteps
    NMeasurements= sum(Sim.SlabHistogramSeries[0,start:stop,1].!=0.0)
    AvgHist[:,:,end] .= OffsetArrays.no_offset_view((reduce(+, Sim.SlabHistogramSeries[xaxis,start:stop,:], dims=2)./(NMeasurements))[:,1,:])

    # ,  title="slab density of $(Sim.SimulationName)"
    fig = Plots.plot(dpi=300,  ylabel= "avg. density\n"* "  [kg/L]" , xlabel= "z-Axis [Å]"  ,  title="$(Sim.SimulationName)",palette=:darktest, seriescolor=:darktest, color=:darktest, markercolor=:darktest, linecolor=:darktest, zcolor=:darktest, xticks=xticks, colorbar_title="window No.", levels=tmp) #,size=(600,300), figsize=(12cm, 4.5cm) Plots.palette(:darktest, 10)
    for i in 1:tmp
        if (i)*Windowlength < Sim.EquilibrationTime continue end
        Plots.plot!(axes(AvgHist)[1][2:end-1], AvgHist[:,1,i][2:end-1],  label="",  line_z=i, marker_z=i, levels=tmp, c=:darktest)
    end

    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_AvgSlabHist_Evolution.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_AvgSlabHist_Evolution.pdf")
    return fig
end

function plotDensityHistogram(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

    x = collect(0:1999)
    fig = Plots.plot(x, Sim.DensityHist, xlim=(10,1000), ylim=(0,maximum(Sim.DensityHist[10:end])), xlabel="m [kg/L]")
    println("In Plotting.")
    println(Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_DensityHist.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_DensityHist.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_DensityHist.pdf")
end