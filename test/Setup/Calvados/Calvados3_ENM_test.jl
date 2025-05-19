using JSON

### test if all bonds are valid; not tested if all bonds that should be there are there....

@testset "Calvados3 ENM" begin 

BasePath = "$SetupTestPath/ENM_Test/"

rm(BasePath; force=true, recursive=true)
mkpath(BasePath)


DomainDict= Dict("RS31" => [[1,70], [90,150]], "RS31a" => [[1,75], [90,140]]) 
Proteins = ["RS31", "RS31", "RS31a"]
ProteinJSON= Dict("RS31" =>"/localscratch/test/fold_rs31/fold_rs31_full_data_0.json","RS31a" =>"/localscratch/test/Alpha_Fold_datas/fold_rs31a/fold_rs31a_full_data_0.json" )
ProteinToCif= Dict("RS31" =>"/localscratch/test/fold_rs31/fold_rs31_model_0.cif","RS31a" =>"/localscratch/test/Alpha_Fold_datas/fold_rs31a/fold_rs31a_model_0.cif" )


HPSAnalysis.RewriteCifToPDB(BasePath,ProteinToCif, Proteins )



for rcut in [7.0,8.0,10.0]
    ConstraintDict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = rcut, plDDTcut=90.0)
    RS31_r = getindex.(ConstraintDict["RS31"],3)
    RS31a_r = getindex.(ConstraintDict["RS31a"],3)
    @test all(RS31_r .<=rcut)
    @test all(RS31a_r .<=rcut) 
end


ConstraintDict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0, pae_cut=1.85)
RS31_i = getindex.(ConstraintDict["RS31"],1)
RS31_j = getindex.(ConstraintDict["RS31"],2)
dom1 = DomainDict["RS31"][1]
dom2 = DomainDict["RS31"][2]
 
isinbounds(x, dom1, dom2) = (x>=dom1[1]&& x<=dom1[2] )|| (x>=dom2[1]&& x<=dom2[2])

#check if ids are within bounds
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31_i ))
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31_j ))


RS31a_i = getindex.(ConstraintDict["RS31a"],1)
RS31a_j = getindex.(ConstraintDict["RS31a"],2)
dom1 = DomainDict["RS31a"][1]
dom2 = DomainDict["RS31a"][2]

#check if ids are within bounds
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31a_i ))
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31a_j ))


RS31_pae  =JSON.parsefile(ProteinJSON["RS31"])["pae"]
RS31a_pae =JSON.parsefile(ProteinJSON["RS31a"])["pae"]


### check if pae values are valid
check(i,j,_, pae, cut) =pae[i][j]<cut
@test reduce(*, map(x->check(x..., RS31_pae, 1.85),  ConstraintDict["RS31"]))
@test reduce(*, map(x->check(x..., RS31a_pae, 1.85),  ConstraintDict["RS31a"]))

CifData = Dict()
RS31_plDDT  = [parse.(Float64,line[15]) for line in split.(strip.(readlines(ProteinToCif["RS31"] ))) if line[1]=="ATOM" && line[4]=="CA"]
RS31a_plDDT = [parse.(Float64,line[15]) for line in split.(strip.(readlines(ProteinToCif["RS31a"]))) if line[1]=="ATOM" && line[4]=="CA"]
for (plDDTcut, pae_cut) in zip([80.0,90.0,92.0], [1.7, 1.85,2.0])
    local ConstraintDict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=plDDTcut, pae_cut=pae_cut)

    @test reduce(*, map(x->check(x..., RS31_pae,pae_cut),  ConstraintDict["RS31"]))
    @test reduce(*, map(x->check(x..., RS31a_pae,pae_cut),  ConstraintDict["RS31a"]))

    local RS31_i  = getindex.(ConstraintDict["RS31"],1)
    local RS31_j  = getindex.(ConstraintDict["RS31"],2)
    local RS31a_i = getindex.(ConstraintDict["RS31a"],1)
    local RS31a_j = getindex.(ConstraintDict["RS31a"],2)
    
    @test reduce(*, map(i-> RS31_plDDT[i]>plDDTcut, RS31_i))
    @test reduce(*, map(i-> RS31_plDDT[i]>plDDTcut, RS31_j))
    @test reduce(*, map(i-> RS31a_plDDT[i]>plDDTcut, RS31a_i))
    @test reduce(*, map(i-> RS31a_plDDT[i]>plDDTcut, RS31a_j))
end


###@TODO: Add manual test for full 
##DomainDict= Dict("RS31" => [[1,20], [90,100]], "RS31a" => [[1,30], [90,105]]) 
end