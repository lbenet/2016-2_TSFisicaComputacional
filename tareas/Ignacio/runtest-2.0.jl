# Este archivo incluye los tests del módulo AD
include("AutomDiff-2.0.jl")
using Base.test
using AD


@test NumDual(5.,5).derivada == 5.
@test NumDual(3).funcion == 3 && NumDual(3).derivada == 0
@test NumDual(3.).funcion == 3. && NumDual(3.).derivada == 0


@test NumDual(2,2)+2 == NumDual(4,2)
@test 2 + NumDual(2,2) == NumDual(4,2)
@test 2*NumDual(2,0) == NumDual(2,0)*2
@test 1/NumDual(2,0) == NumDual(1/2,0)
@test NumDual(1,0)^2 == NumDual(1,0)


