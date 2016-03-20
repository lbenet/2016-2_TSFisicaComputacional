
#Aquí van los tests que se implementaron
include("AutomDiff.jl")

using Base.Test

using AD

# Aqui se incluyen las pruebas necesarias
Dual(3.0,3)-Dual(2,0.5)

Dual(1,3)*Dual(2,5)

Dual(2,1)/Dual(4,2)

Dual(2,3)^2

Dual(3,-4)

Dual(3,2)+1

1+Dual(3,2)

f_test(x) = 3*x^3 - 2

f_test(xdual(1))

gtest(x) = (3x^2-8x+5)/(7x^3-1)

gtest( xdual(1) )


