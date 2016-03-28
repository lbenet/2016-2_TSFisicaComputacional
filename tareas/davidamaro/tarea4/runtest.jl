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
