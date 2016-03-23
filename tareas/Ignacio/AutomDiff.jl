# El modulo en cuestion implementa diferenciacion automatica en Julia,
# definiendo a los duales y a sus operaciones
# Autores: Yuriko Yamamoto, Ignacio Vargas
# Fecha: 22 de marzo, 2016

# La siguiente instrucción sirve para *precompilar* el módulo
__precompile__(true)

module AD
    import Base: +, -, *, /, ^

    export NumDual, xdual

    immutable NumDual{T <: Real}
        funcion::T
        derivada::T
    end

    import Base.convert
    import Base.promote_rule
    import Base.+

    NumDual(a,b) = NumDual(promote(a,b)...)
    NumDual(a) = NumDual(a, zero(a))

    convert{T<:Real}(::Type{NumDual{T}}, a::T) = Dual(a)
    convert{T<:Real, S<:Real}(::Type{NumDual{T}}, a::S) = NumDual(convert(T,a))

    promote_rule{T<:Real, S<:Real}(::Type{NumDual{T}}, ::Type{S}) =
        NumDual{(promote_type)(T,S)}

    function xdual(x0)
        NumDual(x0,1)
    end

    +(a::NumDual, b::NumDual) = NumDual(a.funcion + b.funcion, a.derivada + b.derivada)
    +(a::Real, b::NumDual) = NumDual(a + b.funcion, b.derivada)
    +(a::NumDual, b::Real) = NumDual(a.funcion + b, a.derivada)
    -(a::NumDual, b::NumDual) = NumDual(a.funcion - b.funcion, a.derivada - b.derivada)
    -(a::Real, b::NumDual) = NumDual(a - b.funcion, b.derivada)
    -(a::NumDual, b::Real) = NumDual(a.funcion - b, a.derivada)
    *(a::NumDual, b::NumDual) = NumDual(a.funcion * b.funcion, a.funcion*b.derivada + b.funcion*a.derivada)
    *(a::Real, b::NumDual) = NumDual(a * b.funcion, a*b.derivada)
    *(a::NumDual, b::Real) = NumDual(a.funcion * b, b*a.derivada)
    /(a::NumDual, b::NumDual) = NumDual(a.funcion / b.funcion, (a.derivada - ((a.funcion/b.funcion)*b.derivada)) / b.funcion)
    /(a::Real, b::NumDual) = NumDual(a / b.funcion, ((a/b.funcion)*b.derivada) / b.funcion)
    /(a::NumDual, b::Real) = NumDual(a.funcion / b, a.derivada / b)
    ^(a::NumDual, b::Integer) = NumDual(a.funcion^b, (b*(a.funcion^(b-1))) * a.derivada)
    ^(a::NumDual, b::Real) = NumDual(a.funcion^b, (b*(a.funcion^(b-1))) * a.derivada)
    ==(a::NumDual, b::NumDual) = a.funcion==b.funcion && a.derivada==b.derivada

end
