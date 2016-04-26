include("TaylorSeries.jl")
using Base.Test
using TS

A=(0,1,7,0.0,1.0,7.0,pi,-1,-7)
for (a in A)
    for(b in A)
        d=Taylor([a,b]);
        c=a
        if(a!=0)
            @test (d+d-d*d/d)==d
            @test (d+a-a*a/a)==d
            @test (c+d-d*d/d)==Taylor(c)
        else
            @test ((d+d-d)*d)==d*d
            @test ((d+a-a)*a)==d*a
            @test ((c+d-d)*d)==Taylor(c*d)
        end
    end
end