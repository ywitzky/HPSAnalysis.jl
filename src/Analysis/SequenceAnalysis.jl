using ..BioData

### taken from https://pubs.aip.org/aip/jcp/article/148/12/123305/1023768/Sequence-charge-decoration-dictates-coil-globule
function SequenceChargeDecorator(charges::Vector{R}) where {R<:Real} 
    N = length(charges)
    SCD = zero(R)
    for m in 2:N
        for n in 1:m-1
            SCD += charges[m]*charges[n]*sqrt(m-n)
        end
    end
    return SCD/N
end

### assumes Calvados2 and pH 7 , H needs to be modified for other pH values
function ConvertSequenceToCharges(Sequence::String, Type=Float32, OneToCharge=BioData.OneToCalvados2Charge, pH=7.0)
    #OneToCharge['H']=1. / ( 1 + 10^(pH-6) )

    charges=zeros(Type, length(Sequence))
    for (i,residue) in enumerate(Sequence)
        charges[i]=OneToCharge[residue]
    end
    return charges
end