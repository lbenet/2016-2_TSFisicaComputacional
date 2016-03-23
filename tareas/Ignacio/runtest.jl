# Este archivo incluye los tests del módulo AD
include("AutomDiff.jl")
using Base.test
using AD

@test Dual(5.,5).der == 5.
@test Dual(3).fun == 3 && Dual(3).der == 0
@test Dual(3.).fun == 3. && Dual(3.).der == 0
d2 = Dual(7)
@test d2.fun == 7 && d2.der == 0
#Correcto. No marca error

@test Dual(2,2)+2 == Dual(4,2)
@test 2 + Dual(2,2) == Dual(4,2)
@test 2*Dual(2,0) == Dual(2,0)*2
@test 1/Dual(2,0) == Dual(1/2,0)
@test Dual(1,0)^2 == Dual(1,0)
d1 = Dual(2.0,3.0)
d2 = Dual(2.1,3.0)

@test d1 + d2 == Dual(4.1, 6)
@test d1 + 2.0 == Dual(4.0, 3.0)
@test 2.0 + d1 == Dual(4.0, 3.0)
pres = Dual(1000, 1000)
res = d1 - d2
@test convertir(redondear((res*pres)))/pres == Dual(-0.1, 0)
@test d1*d2 == Dual(4.2, 12.3)
@test d1^2 == Dual(4, 12)
