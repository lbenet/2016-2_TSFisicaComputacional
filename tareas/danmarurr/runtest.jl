include("AutomDiff.jl")
using Base.Test
using AD


a = Dual(1)
b = Dual(1.0,0.0)

@test a.der == 0
@test xdual(1.2).der == 1
@test xdual(1.2).fun == 1.2
@test xdual(3).der == 1.0
@test xdual(3.0).fun == 3

a = Dual(1)
b = xdual(2)

@test a+b == Dual(3, 1)
@test a-b == Dual(-1, -1)
@test a*b == Dual(2, 1)
@test a/b == Dual(1//2, -1//4)
@test b/a == Dual(2, 1)
@test 2*a == Dual(2.0,0.0)
@test a + 1 == Dual(2,0)
@test a^3 == Dual(1,0)
@test a/5 == Dual(1/5, 0)
@test 5/b == Dual(5/2, -5/4)
@test 1 + b == Dual(3,1)