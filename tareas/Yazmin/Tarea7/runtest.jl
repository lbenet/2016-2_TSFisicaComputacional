
#Aquí van los tests que se implementaron con operaciones entre Taylor's.

include("AutomDiff.jl")

using Base.Test

using TYLR

#Operaciones con Taylor's

@test Taylor([1,2.0],0) == Taylor([1.0,2.0],1)

@test Taylor([5,7.0],2) == Taylor([5.0,7.0,0.0],2)

@test Taylor([2,2],3) + Taylor([1,3],5) == Taylor([3,5],3)

@test Taylor([6,4],3) - Taylor([4,1],5) == Taylor([2,3],3)

@test Taylor([3,4],3) * Taylor([4,1],4) == Taylor([12.0,19.0,4.0,0.0],3)

@test Taylor([4,6],3) / Taylor([2,3],3) == Taylor([2.0,0.0,0.0],2)

