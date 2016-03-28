__precompile__(true)

module AD
	import Base: +, -, *, /, ^, ==
	
	export Dual, xdual
"""
Dual

Es un tipo parametrizado por subtipos Reales. Posee entradas fun y der.
"""
type Dual{T<:Real}
    # Se declara el parámetro y las entradas fun y der.
    fun::T
    der::T
end


# Con este método se promueven fun y der para ser del mismo tipo usando
# promote y luego usando un splat dentro del tipo Dual
Dual(fun,der) = Dual(promote(fun, der)...)

# Aqui se define un método que garantiza que el dual de un número cumple lo requerido
Dual(a::Real) = Dual(a,0...)

# Aqui se define la función `xdual`, que se usará para identificar la variable independiente
function xdual(x)
    return Dual(x, 1...) # Como consiste en la variable independiente, es la identidad.
end

#Se definen los métodos necesarios.
# Primero la suma entre duales
+(a::Dual, b::Dual) = Dual((a.fun + b.fun),(a.der + b.der)...)
# La suma entre una constante y una dual.
+(a::Real, b::Dual) = Dual(a) + b
+(a::Dual, b::Real) = b + a
# La resta entre duales y duales y constantes
-(a::Dual, b::Dual) = Dual((a.fun - b.fun),(a.der - b.der)...)
-(a::Real, b::Dual) = Dual(a)-b
-(a::Dual, b::Real) = a - Dual(b)
# La multiplicación entre duales y duales y constantes.
*(a::Dual, b::Dual) = Dual((a.fun * b.fun),(a.fun * b.der + b.fun * a.der)...)
*(a::Real, b::Dual) = Dual(a) *b
*(a::Dual, b::Real) = b*a
# La división entre duales.
/(a::Dual, b::Dual) = Dual((a.fun / b.fun),((a.der-(a.fun / b.fun)*b.der)/b.fun)...)
# La exponenciación, para esta, Julia nos pide que definamos primero la exponenciación
# con un exponente entero y luego ya con un número Real.
^(a::Dual,ex::Integer) = Dual((a.fun^ex), (ex*a.fun^(ex-1)*a.der)...)
^(a::Dual,ex::Real) = Dual((a.fun^ex), (ex*a.fun^(ex-1)*a.der)...)
# Finalmente, los operadores unitarios y de comparacion.
+(a::Dual) = a
-(a::Dual) = Dual(-a.fun, -a.der...)
==(a::Dual, b::Dual) = (a.fun == b.fun && a.der == b.der) ? true : false

end
