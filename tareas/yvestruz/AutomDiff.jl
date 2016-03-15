#=
Este módulo contiene las bases para la diferenciación automática en Julia. Permite calcular de manera exacta la derivada de cualquier función racional real, introduciendo un nuevo tipo: los duales.

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

^(a::Dual, b::Int) = Dual(a.fun^b,b*(a.fun^(b-1))*a.der)
^(a::Dual, b::Real) = Dual(a.fun^b,b*(a.fun^(b-1))*a.der)

==(a::Dual, b::Dual) = (a.fun == b.fun) && (a.der == b.der)



end
