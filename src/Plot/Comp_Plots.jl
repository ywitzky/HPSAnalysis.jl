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
