# Common Slab analysis methods


### Computation of slab histograms analysis

Here is a common way how slab simulations are analysed.

```julia
import HPSAnalysis as HPS

prot= HPS.initData("/Path/To/Data/"; Reparse=true, LoadAll=true, HOOMD=true)

prot.EquilibrationTime=...

### compute RG and it's auto correlation 
HPS.computeCOM!(prot)
HPS.computeRGSeries!(prot)
HPS.computeRGCorrelationTime(prot)
HPS.plotRGAutocorr(prot)

HPS.computeClustersByChainCOM(prot)
### computes slab histograms in range 
### Sim.ClusterRange=prot.EquilibrationTime:prot.RGMeasureStep:prot.NSteps
HPS.computeSlabHistogram(prot)
HPS.plotAvgSlabDensity(prot; Windowlength=1000)
HPS.plotAvgSlabDensityEvolution(prot; Windowlength=100)

HPS.SaveAllData(prot)


```