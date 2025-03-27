using HPSAnalysis,Test
using HPSAnalysis.BioData
using GSDFormat

Sequences=["MRPVFV","MRPVF","MRPV","MRP"]
SimulationType="Calvados2"
pH=7.5

@testset "DetermineCalvados2AtomTypes" begin
    AtomTypes_test=Set(join(Sequences))
    AaToId_test = Dict{Char,Int32}()
    for (index, value) in enumerate(AtomTypes_test)
        AaToId_test[value]=index
    end
    Long=Set()
    OneToCharge_test=Dict() 
    OneToMass_test=deepcopy(BioData.AaToWeight)
    OneToSigma_test=Dict() 
    OneToLambda_test=Dict() 
    OneToHPSDihedral0110_test=deepcopy(BioData.OneToHPSDihedral0110)
    OneToHPSDihedral1001_test=deepcopy(BioData.OneToHPSDihedral1001)
    LongAtomTypesToRes_test=Dict{Char,Tuple{Char,Bool}}()

    if SimulationType=="HPS-Alpha"
        OneToCharge_test = deepcopy(BioData.OneToHPSCharge)
        OneToLambda_test = deepcopy(BioData.OneToHPSUrryLambda)
        OneToSigma_test  = deepcopy(BioData.OneToHPSCalvadosSigma)
    elseif SimulationType=="Calvados2"
        AaToId_test['a']=Int32(6)
        AaToId_test['b']=Int32(7)
        AaToId_test['c']=Int32(8)
        AaToId_test['d']=Int32(9)
        Long=Set(['a','b','c','d'])
        LongAtomTypesToRes_test['a']=('M',1)
        LongAtomTypesToRes_test['b']=('V',0)
        LongAtomTypesToRes_test['c']=('F',0)
        LongAtomTypesToRes_test['d']=('P',0)
    else
        OneToCharge_test = deepcopy(BioData.OneToHPSCharge)
        OneToLambda_test = deepcopy(BioData.OneToCalvados2Lambda)
        OneToSigma_test = deepcopy(BioData.OneToHPSCalvadosSigma)
    end

    LongAtomTypes_test=union(AtomTypes_test, Long)

    if SimulationType=="Calvados2"
        OneToCharge_test = deepcopy(BioData.OneToHPSCharge)
        OneToLambda_test = deepcopy(BioData.OneToCalvados2Lambda)
        OneToSigma_test  = deepcopy(BioData.OneToHPSCalvadosSigma)
        OneToCharge_test['H'] = 1. / ( 1 + 10^(pH-6) ) 
        for e in LongAtomTypes_test
            if ~(e in keys(OneToCharge_test))
                (AA,front)=LongAtomTypesToRes_test[e]
                OneToCharge_test[e] = front ? OneToCharge_test[AA] +1 :  OneToCharge_test[AA] -1
                OneToMass_test[e] = front ? OneToMass_test[AA] +2.0 :  OneToMass_test[AA] +16
                OneToSigma_test[e] = OneToSigma_test[AA]
                OneToLambda_test[e] = OneToLambda_test[AA]
                OneToHPSDihedral0110_test[e] = OneToHPSDihedral0110_test[AA]
                OneToHPSDihedral1001_test[e] = OneToHPSDihedral1001_test[AA]
            end
        end
    end

    IdToAa_test=Dict((v=>k) for (k,v) in AaToId_test)

    (AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)=HPSAnalysis.Setup.DetermineCalvados2AtomTypes(Sequences,SimulationType,pH)

    @test (AtomTypes_test==AtomTypes)
    @test (LongAtomTypes_test==LongAtomTypes)
    @test (AaToId_test==AaToId)
    @test (IdToAa_test==IdToAa)
    #@test (ResToLongAtomType_test==ResToLongAtomType)
    @test (LongAtomTypesToRes_test==LongAtomTypesToRes)
    @test (OneToCharge_test==OneToCharge)
    @test (OneToMass_test==OneToMass)
    @test (OneToSigma_test==OneToSigma)
    @test (OneToLambda_test==OneToLambda)
    @test (OneToHPSDihedral0110_test==OneToHPSDihedral0110)
    @test (OneToHPSDihedral1001_test==OneToHPSDihedral1001)
end

@testset "Yukawa Interaction" begin
    SimulationType="Calvados2"
    Temperature=300
    SaltConcentration=0.5
    Temp=300
    ionic=0.5

    if SimulationType=="Calvados2"
        e = 1.6021766### Charge of electron
        e_0 = 8.854188### vacuum permitivity
        NA = 6.022#14086# 1/mol Avogadro constant
        kb = 0.00831446262
        kT = kb*Temp

        epsilon_r=5321.0/Temp+233.76-0.9297*Temp+1.417*10.0^(-3)*Temp^2-8.292*10.0^(-7)*Temp^3
        lamb = e^2/(4.0*pi*e_0*epsilon_r)*NA*1000/kT
        kappa=sqrt((8.0*pi*lamb*ionic*NA/10))

        ϵ_r_test=epsilon_r
        κ_test=kappa
    else
        ϵ_r_test=80.0
        κ_test=10.0
    end

    ϵ_r, κ=HPSAnalysis.Setup.DetermineYukawaInteractions(;SimulationType,Temperature,SaltConcentration)

    @test (ϵ_r_test≈ϵ_r)
    @test isapprox(κ_test,κ; atol=1e-5)
end
