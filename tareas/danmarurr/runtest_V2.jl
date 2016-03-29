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

d1 = xdual(0.5)
d2 = xdual(1.5)

@test sin(d1) == Dual(sin(d1.fun), cos(d1.fun)*d1.der)
@test cos(d1) == Dual(cos(d1.fun), - sin(d1.fun)*d1.der)

Vec_Func = [(:sin, :cos), (:cos, :(x -> -sin(x))), (:tan, :(x -> (sec(x))^2)), (:cot, :(x -> -(csc(x))^2)), 
    (:sec, :(x -> sec(x)*tan(x))), (:csc, :(x -> -csc(x)*cot(x))), (:sinh, :cosh), (:cosh, :sinh), 
    (:tanh, :(x -> (sech(x))^2)), (:coth, :(x -> -(csch(x))^2)), (:asin, :(x -> 1/sqrt(1-x^2))), 
    (:acos, :(x -> -1/sqrt(1-x^2))), (:atan, :(x -> 1/(1+x^2))), (:acot, :(x -> -1/(1+x^2))),
    (:asec, :(x -> 1/(sqrt(x^2-1)*x))), (:acsc, :(x -> -1/(sqrt(-1+x^2)*x))), (:asinh, :(x -> 1/sqrt(1+x^2))),
    (:acosh, :(x -> 1/sqrt(x^2-1))), (:atanh, :(x -> 1/(1-x^2))), (:acoth, :(x -> 1/(1-x^2))),
    (:asech, :(x -> 1/(x*sqrt(1-x^2)))), (:acsch, :(x -> -1/(x*sqrt(1+x^2)))),
    (:sqrt, :(x -> 1/(2*sqrt(x)))), (:exp, :exp), (:cbrt, :(x -> 1/(3*x^(2/3))))]

#Creamos un arreglo para generar las posiciones de las funciones que aceptan valores menores que 1
arr1 = Int64[]
for r in 1:25
    (r == 15) | (r == 18) | (r==16) | (r == 20)? nothing : push!(arr1, r)
end
#Creamos un arreglo para generar las posiciones de las funciones que aceptan valores mayores que 1
arr2 = Int64[]
for r in 1:25
    (r == 15) | (r == 18) | (r==16) | (r == 20)? push!(arr2, r) : nothing 
end


function probar_en_dual(d::Dual,arr::Array{Int64,1})  
    for r in arr
        fun = Vec_Func[r][1]
        der = Vec_Func[r][2]

        e = quote
            p = ($fun)($d) == Dual(($fun)($d.fun), (($der)($d.fun))*$d.der) 
        end
        #@show r
        @eval $e
        @test p
        #println(p)
    end            
end

probar_en_dual(d1, arr1)
probar_en_dual(d2, arr2)
     
