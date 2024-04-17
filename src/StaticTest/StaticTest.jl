using StaticArrays, BenchmarkTools

a= [i for i in 1:100]
b= MMatrix{10,10,Int64}(reshape(@view a,(10,10)))

function test(a)
    r = zeros(10)
    for (i,(b,c)) in enumerate(zip(1:10:90, 9:10:99))
        r[i] = sum(a[b:c])
    end
    return r
end

function test(a)
    r = zeros(10)
    for (i,(b,c)) in enumerate(zip(1:10:90, 9:10:99))
        r[i] = sum(a[b:c])
    end
    return r
end


@btime test(a)
@btime test2(b)


mutable struct test13{N, B}
    equichain::B
    a::MMatrix{N,N, Int32}
end


a = test13{3}(MMatrix{3,3, Int32}(0,0,0,0,0,0,0,0,0))

a.a = MMatrix{3,3, Int32}(0,0,0,0,0,0,0,0,0)
