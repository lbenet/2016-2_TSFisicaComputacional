# El modulo en cuestion implementa diferenciacion automatica en Julia,
# definiendo a los duales y a sus operaciones
# Autores: Yuriko Yamamoto, Ignacio Vargas
# Fecha: 29 de marzo, 2016

# La siguiente instrucción sirve para *precompilar* el módulo
__precompile__(true)

module AD
    import Base: +, -, *, /, ^

    export Dual, xdual

    type Dual{T <: Real}
        fun::T
        der::T
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
    
    ^(a::Dual, b::Integer) = Dual(a.fun^b, (b*(a.fun^(b-1))) * a.der)
    ^(a::Dual, b::Real) = Dual(a.fun^b, (b*(a.fun^(b-1))) * a.der)
    
    ==(a::Dual, b::Dual) = a.fun==b.fun && a.der==b.der
    
    
end
