#=
Este módulo contiene las bases para la diferenciación automática en Julia. Permite calcular de manera exacta la derivada de cualquier función racional real y para ciertas funciones especiales, introduciendo un nuevo tipo: los duales.

Autores:
Luis Yves Villegas Aguilar
David Leonardo Galicia Praskauer
Arturo Cid Prati

15 de Marzo de 2016
=#


__precompile__(true)

module AD
    import Base: +, -, *, /, ^

    export Dual, xdual

    
type Dual{T<:Real}
    # código: 
    fun :: T
    der :: T
end

Dual(a,b) = Dual(promote(a,b)...)
Dual(a::Real) = Dual(a,0)

function xdual(x0)
        return Dual(x0,1)
end

import Base: +, -, *, /, ^, ==

# Aqui se implementan los métodos necesarios para cada función
+(a::Dual, b::Dual) = Dual(a.fun + b.fun, a.der + b.der)
+(a::Dual, b::Real) = Dual(a.fun + b, a.der)
+(b::Real, a::Dual) = Dual(a.fun + b, a.der)
+(a::Dual) = Dual(a.fun,b.der)

-(a::Dual, b::Dual) = Dual(a.fun - b.fun, a.der - b.der)
-(a::Dual, b::Real) = Dual(a.fun - b, a.der)
-(b::Real, a::Dual) = Dual(b - a.fun, a.der)
-(a::Dual) = Dual(-1*a.fun,-1*a.der)

*(a::Dual, b::Dual) = Dual(a.fun * b.fun, (a.fun * b.der) + (b.fun * a.der))
*(a::Real, b::Dual) = Dual(a*b.fun, a*b.der)
*(b::Dual, a::Real) = Dual(a*b.fun, a*b.der)

/(a::Dual, b::Dual) = Dual(a.fun/b.fun, (a.der - (a.fun/b.fun)*b.der)/b.fun)
/(b::Dual, a::Real) = Dual((b.fun)/a,(b.der)/a)
/(a::Real, b::Dual) = Dual(a/(b.fun),a/(b.der))  #ésta fue la que añadí

^(a::Dual, b::Integer) = Dual(a.fun^b,b*(a.fun^(b-1))*a.der)
^(a::Dual, b::Real) = Dual(a.fun^b,b*(a.fun^(b-1))*a.der)

==(a::Dual, b::Dual) = (a.fun == b.fun) && (a.der == b.der)

import Base: sqrt, cbrt

sqrt(a::Dual) = Dual(sqrt(a.fun),(1//2)*(a.der/sqrt(a.fun)))
cbrt(a::Dual) = Dual(cbrt(a.fun),(1//3)* (a.der/(cbrt(a.fun)^2)))

import Base: exp, log

exp(a::Dual) = Dual(exp(a.fun), a.der*exp(a.fun))
log(a::Dual) = Dual(log(a.fun), a.der/a.fun)

#funciones trigonométricas con argumentos en radianes

import Base: sin, cos, tan, cot, sec, csc

sin(a::Dual) = Dual(sin(a.fun), a.der*cos(a.fun))
cos(a::Dual) = Dual(cos(a.fun),-a.der*sin(a.fun))

tan(a::Dual) = Dual(tan(a.fun), a.der*sec(a.fun)*sec(a.fun))
cot(a::Dual) = Dual(cot(a.fun), -1*a.der*csc(a.fun)*csc(a.fun))

sec(a::Dual) = Dual(sec(a.fun), a.der*tan(a.fun)*sec(a.fun))
csc(a::Dual) = Dual(csc(a.fun), a.der*(-1*cot(a.fun))*csc(a.fun))

#funciones trigonométricas hiperbólicas

import Base: sinh, cosh, tanh, coth, sech, coth, csch

sinh(a::Dual) = Dual(sinh(a.fun), a.der*cosh(a.fun))
cosh(a::Dual) = Dual(cosh(a.fun), a.der*sinh(a.fun))

tanh(a::Dual) = Dual(tanh(a.fun), a.der*sech(a.fun)*sech(a.fun))
coth(a::Dual) = Dual(coth(a.fun), -1*a.der*csch(a.fun)*csch(a.fun))

sech(a::Dual) = Dual(sech(a.fun), a.der*tanh(a.fun)*(-1*sech(a.fun)))
csch(a::Dual) = Dual(csch(a.fun), a.der*(-1*coth(a.fun))*csch(a.fun))

end
