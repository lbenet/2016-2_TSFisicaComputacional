# Este archivo incluye los tests del módulo AD
include("AutomDiff-t5.jl")
using Base.test
using AD


@test xdual(2)==AD.Dual{Int64}(2,1)
@test f(xdual(2))==AD.Dual{Float64}(0.9934213368955197,1)
@test newton1D(W6, 2.2)==2.0
@test newton1D(W6, 2.45)==6.0