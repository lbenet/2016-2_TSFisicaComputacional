# El modulo en cuestion implementa diferenciacion automatica en Julia,
# definiendo a los duales y a sus operaciones
# Autores: Yuriko Yamamoto, Ignacio Vargas
# Fecha: 22 de marzo, 2016

# La siguiente instrucción sirve para *precompilar* el módulo
__precompile__(true)

module AD
    import Base: +, -, *, /, ^

    export Dual, xdual

    """Definición de los duales, donde 'fun' es la variable de la funcion
    evaluada en x_0, y 'der' es la variable de la derivada en x_0
    ...
    """
    type Dual{T<:Real} #<: Real
        fun :: T
        der :: T
    end

    import Base.convert
    import Base.promote_rule
    import Base.+

    Dual(a,b) = Dual(promote(a,b)...)
    Dual(a) = Dual(a, zero(a))

    convert{T<:Real}(::Type{Dual{T}}, a::T) = Dual(a)
    convert{T<:Real, S<:Real}(::Type{Dual{T}}, a::S) = Dual(convert(T,a))

    promote_rule{T<:Real, S<:Real}(::Type{Dual{T}}, ::Type{S}) =
        Dual{(promote_type)(T,S)}
        
    function xdual(x0)
        Dual(x0,1)
    end

    +(a::Dual, b::Dual) = Dual(a.fun + b.fun, a.der + b.der)
    +(a::Real, b::Dual) = Dual(a + b.fun, b.der)
    +(a::Dual, b::Real) = Dual(a.fun + b, a.der)

    -(a::Dual, b::Dual) = Dual(a.fun - b.fun, a.der - b.der)
    -(a::Real, b::Dual) = Dual(a - b.fun, b.der)
    -(a::Dual, b::Real) = Dual(a.fun - b, a.der)

    *(a::Dual, b::Dual) = Dual(a.fun * b.fun, a.fun*b.der + b.fun*a.der)
    *(a::Real, b::Dual) = Dual(a * b.fun, a*b.der)
    *(a::Dual, b::Real) = Dual(a.fun * b, b*a.der)

    /(a::Dual, b::Dual) = Dual(a.fun / b.fun, (a.der - ((a.fun/b.fun)*b.der)) / b.fun)
    /(a::Real, b::Dual) = Dual(a / b.fun, ((a/b.fun)*b.der) / b.fun)
    /(a::Dual, b::Real) = Dual(a.fun / b, a.der / b)

    ^(a::Dual, c::Integer) = Dual(a.fun^c, (c*(a.fun^(c-1))) * a.der)
    ^(a::Dual, c::Real) = Dual(a.fun^c, (c*(a.fun^(c-1))) * a.der)

    redondear(a::Dual) = Dual(round(a.fun), round(a.der))

    convertir(a::Dual) = Dual(Float64(a.fun), round(a.der))

    ==(a::Dual, b::Dual) = a.fun==b.fun && a.der==b.der

end
