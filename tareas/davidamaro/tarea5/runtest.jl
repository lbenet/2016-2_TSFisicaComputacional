# Este archivo incluye los tests del módulo AD
include("AutomDiff.jl")
using Base.Test
using AD

# Potencia
@test Dual(1,0)^2 == Dual(1,0)
@test Dual(2,0)^2 == Dual(4,0)
# Suma
@test Dual(1,2)+2 == Dual(3,2)
@test 2 + Dual(1,2) == Dual(3,2)
# Producto y división
@test 2*Dual(4,0) == Dual(4,0)*2
@test 1/Dual(3,0) == Dual(1/3,0)

@test log(xdual(1)) == Dual(log(1),1/1) # d/dx log(x) en x=1
@test exp(4xdual(4)) == Dual(exp(4*4),4*exp(16)) # d/dx e^4x en x = 4
@test exp(xdual(1)) == Dual(exp(1),exp(1)) # d/dx e^x en x = 1
@test sin(4xdual(2)) == Dual(sin(4*2),4*cos(4*2)) # d/dx sin(4x) en x=2
@test tan(4xdual(4.2)) == Dual(tan(4*4.2),4*sec(4*4.2)^2) # d/dx tan(4x) en x = 4.2
@test sinh(2xdual(1.5)) == Dual(sinh(2*1.5),2*cosh(2*1.5)) # d/dx sinh(2x) en x = 1.5
