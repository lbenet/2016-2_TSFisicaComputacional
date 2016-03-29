include("AutomDiff.jl")
using Base.test
using AD

A=(0,1,7,0.0,1.0,7.0,pi,-1,-7)
for (a in A)
    for(b in A)
        d=Dual(a,b);
        c=xdual(a)
        if(a!=0)
            @test (d+d-d*d/d)^1==d
            @test (d+a-a*a/a)^1==d
            @test (c+d-d*d/d)^1==c
        else
            @test ((d+d-d)*d)^1==d*d
            @test ((d+a-a)*a)^1==d*a
            @test ((c+d-d)*d)^1==c*d
        end
    end
end