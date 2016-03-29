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

@test b^2 == Dual(4, 4)
@test b^(0.5) == Dual(sqrt(2), 0.5/(sqrt(2)))
@test b^(1//2) == Dual(sqrt(2), 0.5/(sqrt(2)))
@test b^pi == Dual(b.fun^pi, pi*b.fun^(pi - 1.0)*b.der)

@test a^2 == Dual(1, 0)
@test a^(0.5) == Dual(sqrt(1), 0)
@test a^(1//2) == Dual(sqrt(1), 0)

#Funcion que prueba en un dual con ayuda de metaprogramming

d1 = xdual(π/2)
function probar_en_dual(d::Dual)  
    for exp in Vec_Func
        fun = exp[1]
        der = exp[2]

        e = quote
            p = ($fun)($d) == Dual(($fun)($d.fun), (($der)($d.fun))*$d.der) 
        end
        
        @eval $e
        #println(p)
        @test p
    end            
end