using SerieTaylor
using Base.Test

@test Taylor([0,5]) + 5 == Taylor([5,5])
@test Taylor([6,6,6]) == Taylor([6,6,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
@test Taylor([7,0]) + Taylor([0,7]) == Taylor([7,7])
@test Taylor(ones(8))/Taylor([1,1]) * Taylor([1,1]) == Taylor(ones(8))
@test Taylor([1.,2.,1.])/Taylor([1.,1.]) == Taylor([1.,1.])

x = Taylor([0,1],15)
z = zero(x)
u = 1.0*one(x)

@test eltype(promote(z,Taylor(u))[1]) == Float64
@test eltype(auxforma(z,u)[1]) == Float64
@test length(auxforma(Taylor(0,5),z)[1]) == 15
@test auxnocero(x) == 1
@test auxnocero(z) == z.orden+1

@test u == 1
@test 0.0 == z
@test x.coeficientes[2] == 1
@test z+1 == u
@test x+x == 2x
@test x-x == z

x2 = Taylor([0,0,1],15)
@test x*x == x2
@test (-x)^2 == x2
@test 1-x2 == (1+x)-x*(1+x)
@test (1-x2)^2 == (1+x)^2 * (1-x)^2
@test log((1-x)^2) == 2*log(1-x)

@test sin(Taylor([0])) == 0
@test cos(Taylor([0])) == 1
@test square(Taylor([2])) == 4
@test ^(Taylor([2]), 2) == 4
@test ^(Taylor([2]), 3) == 8
@test ^(Taylor([2]), 4) == 16

