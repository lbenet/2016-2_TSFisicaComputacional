#Modulo de Diferenciacion Automatica
#El siguiente modulo define la estructura de los objetos "duales" asi como las operaciones entre ellos.
#Autores: Daniel (https://github.com/danmarurr) y Fernanda (https://github.com/FernandaPerez)


__precompile__(true)

module AD
    import Base: +, -, *, /, ^, ==
    
    export Dual, xdual



    type Dual{T<:Real}
    # código: 
    fun :: T
    der :: T
    end
    
    Dual(a, b) = Dual(promote(a, b) ...)
    # Aqui se define un método que garantiza que el dual de un número cumple lo requerido
    Dual(a) = Dual(a, zero(0))
    # Aqui se define la función `xdual`, que se usará para identificar la variable independiente

    function xdual(x0)
        Dual(x0, one(x0))
    end

    # Definiendo operaciones cuando los argumentos son Duales
    +(a::Dual, b::Dual) = Dual(a.fun + b.fun, a.der + b.der)
    -(a::Dual, b::Dual) = Dual(a.fun - b.fun, a.der - b.der)
    *(a::Dual, b::Dual) = Dual(a.fun * b.fun, a.der*b.fun + b.der*a.fun)
    /(a::Dual, b::Dual) = Dual(a.fun / b.fun, (a.der - b.der*(a.fun/b.fun))/b.fun)
    ^{T<:Int}(a::Dual, α::T) = Dual(a.fun^α, α*a.fun^(α - 1)*a.der)
    ^{T<:Float64}(a::Dual, α::T) = Dual(a.fun^α, α*a.fun^(α - 1)*a.der)
    ^{T<:Rational}(a::Dual, α::T) = Dual(a.fun^α, α*a.fun^(α - 1)*a.der)
    ^{T<:Irrational}(a::Dual, α::T) = Dual(a.fun^α, α*a.fun^(α - 1)*a.der)
    ==(a::Dual, b::Dual) = (a.fun == b.fun && a.der == b.der)


    # Definiendo operaciones cuando el argumento de la derecha es un número
    *{T<:Real}(β :: T, a :: Dual) = Dual(a.fun*β, a.der*β    )
    +{T<:Real}(a::Dual, γ::T)= +(a::Dual, Dual(γ))
    -{T<:Real}(a::Dual, γ::T)= -(a::Dual, Dual(γ))
    *{T<:Real}(a::Dual, γ::T)= *(a::Dual, Dual(γ))
    /{T<:Real}(a::Dual, γ::T)= /(a::Dual, Dual(γ))

    # Definiendo operaciones cuando el argumento de la izquierda es un número

    +{T<:Real}(γ::T, a::Dual)= +(Dual(γ), a)
    -{T<:Real}(γ::T, a::Dual)= -(Dual(γ), a)
    *{T<:Real}(γ::T, a::Dual)= *(Dual(γ), a)
    /{T<:Real}(γ::T, a::Dual)= /(Dual(γ), a)


#Importamos todas las funciones para las cuales queremos definir su operacion 
import Base: ^, exp, sqrt, cbrt, sin, cos, tan, cot, sec, csc, sinh, cosh, tanh, coth, sech, csch,
asin,   acos,   atan,   acot,   asec,   acsc,
asinh,  acosh,  atanh,  acoth,  asech,  acsch

#Vector con todos los símbolos asociados a las funciones trigométricas, exponenciales, etc. y sus respectivas
#derivadas
Vec_Func = [(:sin, :cos), (:cos, :(x -> -sin(x))), (:tan, :(x -> (sec(x))^2)), (:cot, :(x -> -(csc(x))^2)), 
    (:sec, :(x -> sec(x)*tan(x))), (:csc, :(x -> -csc(x)*cot(x))), (:sinh, :cosh), (:cosh, :sinh), 
    (:tanh, :(x -> (sech(x))^2)), (:coth, :(x -> -(csch(x))^2)), (:asin, :(x -> 1/sqrt(1-x^2))), 
    (:acos, :(x -> -1/sqrt(1-x^2))), (:atan, :(x -> 1/(1+x^2))), (:acot, :(x -> -1/(1+x^2))),
    (:asec, :(x -> 1/(sqrt(1-x^-2)*x^2))), (:acsc, :(x -> -1/(sqrt(1-x^-2)*x^2))), (:asinh, :(x -> 1/sqrt(1+x^2))),
    (:acosh, :(x -> 1/sqrt(x^2-1))), (:atanh, :(x -> 1/(1-x^2))), (:acoth, :(x -> 1/(1-x^2))),
    (:asech, :(x -> 1/(x*sqrt(1-x^2)))), (:acsch, :(x -> 1/(x*sqrt(1-x^2)))), (:acsch, :(x -> -1/(x*sqrt(1+x^2)))),
    (:sqrt, :(x -> 1/(2*sqrt(x)))), (:exp, :exp), (:cbrt, :(x -> 1/(3*x^(2/3))))]

#Casos especiales: los logarítmos y a^x
log(a::Dual) = Dual(log(a.fun), a.der/a.fun)
log{T<:Real}(b::T, a::Dual) = Dual(log(b,a.fun), a.der/(log(b)*a.fun))
^{T<:Real}(b::T, a::Dual) = Dual(b^a.fun, a.der*log(b)*b^a.fun)

#Loop para crear los nuevos métodos a partir de los símbolos
for r in 1:length(Vec_Func)
    fn = Vec_Func[r][1] #El primer símbolo está asociado a la función
    der = Vec_Func[r][2] #El segundo símbolo está asociado a la derivada de la función
    ex = quote #Creamos una nueva expresión
        function ($fn)(a::Dual) #Definimos fn(a::Dual)
            fun = ($fn)(a.fun)
            derv = ($der)(a.fun)
            return Dual(fun, derv*a.der) #Aplicamos la regla de la cadena
        end
    end
    @eval $ex #Evaluamos la expresión para crear el método
end




end

