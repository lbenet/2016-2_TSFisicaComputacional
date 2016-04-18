__precompile__(true)

module AD
	import Base: +, -, *, /, ^, ==
	import Base: exp, log, sin, cos, tan, sinh, cosh, tanh, asin, 
            acos, atan, asinh, acosh, atanh
	
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
/(a::Dual, b::Real) = a/Dual(b)
/(a::Real, b::Dual) = Dual(a)/b
# La exponenciación, para esta, Julia nos pide que definamos primero la exponenciación
# con un exponente entero y luego ya con un número Real.
^(a::Dual,ex::Integer) = Dual((a.fun^ex), (ex*a.fun^(ex-1)*a.der)...)
^(a::Dual,ex::Real) = Dual((a.fun^ex), (ex*a.fun^(ex-1)*a.der)...)
# Finalmente, los operadores unitarios y de comparacion.
+(a::Dual) = a
-(a::Dual) = Dual(-a.fun, -a.der...)
==(a::Dual, b::Dual) = (a.fun == b.fun && a.der == b.der) ? true : false

# AHORA SE EMPIEZA LA SEGUNDA PARTE DE LA AUTODIFERENCIACIÓN.

# Esta lista de tuplas consiste de las funciones y la derivada correspondiente.
oper = [(:exp, :exp), (:log, :inv), (:sin, :cos), (:cos, :(x -> -sin(x))),
        (:tan, x -> (sec(x))^2), (:sinh, :cosh), (:cosh, :sinh),
        (:tanh, x -> (sinh(x))^2),(:asin, x -> 1/sqrt(1-x^2)), (:acos, x -> -1/sqrt(1-x^2)),
        (:atan, x -> 1/(1+x^2)), (:asinh, x -> 1/sqrt(x^2+1)), 
        (:acosh, x -> 1/(sqrt(x-1)sqrt(x+1))), (:atanh, x -> 1/(1-x^2))];

# En este for se itera sobre la lista anterior y se forman expresiones distintas para cada
# función.
for (f,fd) in oper # f es la función y fd la derivada
    ex = quote
        function  ($f)(a::Dual) # Aquí se define el método f para una dual
            fu = ($f)(a.fun) # Sabemos que la primera entrada de la dual
            # es la misma función.
            de = a.der*($fd)(a.fun) # Para la derivada, tenemos la regla de
            # la cadena que nos dice que mulipliquemos la derivada de la dual
            # que estamos metiendo por la fd correspondiente valuada en a.fun
            return Dual(fu, de)
        end
    end
    @eval $ex # Evaluamos cada una de las expresiones que formamos.
end
end

