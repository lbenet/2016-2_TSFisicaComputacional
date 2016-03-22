
#Aquí van los tests que se implementaron con operaciones entre duales.
include("AutomDiff.jl")

using Base.Test

using AD


@test Dual(3) == Dual(3,0...) 

@test Dual(3.0,3) - Dual(2,0.5) == Dual(3.0-2.0,3-0.5 ...)

@test Dual(1,3)*Dual(2,5) == Dual(2,11...) 

@test Dual(2,1)/Dual(4,2) == Dual(0.5,0.0...)

@test Dual(2,3)^2 == Dual(4,12...)

@test Dual(3,2)+1 == Dual(4,2...)

@test 1+Dual(3,2) == Dual(4,2...)

@test xdual((4)/2) == xdual(2)

@test xdual((10*5)) == xdual(50)

@test xdual((44-4)) == xdual(40)

@test xdual(7) == Dual(7,1...)

f_test(x) = 3*x^3 - 2

@test f_test(xdual(1)) == Dual(1,9...)

gtest(x) = (3x^2-8x+5)/(7x^3-1)

@test gtest( xdual(1) ) == Dual(0.0,-0.3333333333333333...)


