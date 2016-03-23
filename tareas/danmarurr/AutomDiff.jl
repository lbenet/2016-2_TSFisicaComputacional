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


end