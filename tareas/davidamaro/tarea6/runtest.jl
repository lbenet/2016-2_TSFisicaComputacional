# Este archivo incluye los tests del módulo AD
include("Automata.jl")
using Base.Test
using AD

# Multiplicación, división, resta, suma e igualdad.

@test AD.Taylor(0,[1,2,3])+AD.Taylor(0,[2,3,4]) == AD.Taylor(0,[3,5,7])
@test AD.Taylor(0,[1,2,3])-AD.Taylor(0,[2,3,4]) == AD.Taylor(0,[-1,-1,-1])
@test AD.Taylor(0,[0,0,1])/AD.Taylor(0,[0,1]) == AD.Taylor(0,[0,1,0])
@test AD.Taylor(0,[0,0,1,0])/AD.Taylor(0,[0,0,1]) == AD.Taylor(0,[1,0,0,0])
@test AD.Taylor(0,[0,1,1])*AD.Taylor(0,[0,1]) == AD.Taylor(0,[0,0,1,1])

# Funciones elementales

# cos(x) = 1 - x^2/2!+ x^4/4!
cos(AD.Taylor(0,[0,1//1,0,0,0])) == AD.Taylor(0,[1//1,0//1,-1//2,0//1,1//24])
# sin(x) = x -x^3/3!
sin(AD.Taylor(0,[0,1//1,0,0,0])) == AD.Taylor(0,[0//1,1//1,0//1,-1//6,0//1])
# log(x+1) = x-x^2/2+x^3/3-x^4/4+x^5/5
log(AD.Taylor(0,[1//1,1//1,0,0,0,0])) == AD.Taylor(0,[0//1,1//1,-1//2,1//3,-1//4,1//5])
# x^(1/2) en x = 1 es aproximadamente
# 1+(x-1)/2-(x-1)^2/8+(x-1)^3/16-(x-1)^4*5/182
AD.Taylor(1//1,[1//1,1//1,0,0,0])^(1//2) == AD.Taylor(1//1,[1//1,1//2,-1//8,1//16,-5//128])
