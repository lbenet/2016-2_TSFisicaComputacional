# Aqui viene una explicación de lo que se hace en el módulo, los autores y la fecha

# La siguiente instrucción sirve para *precompilar* el módulo
__precompile__(true)

module AD
    import Base: +, -, *, /, ^, ==
    
    export Dual, xdual, cdual
    doc"""
    Definición de los duales, donde
    
    1. `fun` representa el valor de la __función__ evaluada en $x_0$.
    2. `der` el valor de la __derivada__ evaluada en $x_0$.
    
    """
    
    type Dual{T<:Real}
        fun :: T
        der :: T
    end

    # Referencia: La notebook de la clase donde se vió el tema.
    Dual(a, b) = Dual(promote(a,b)...)

    +(a::Dual, b::Dual) = Dual(a.fun + b.fun, a.der + b.der)
    +(a::Dual, b::Real) = a+cdual(b)
    +(b::Real, a::Dual) = a+cdual(b)
    -(a::Dual, b::Dual) = Dual(a.fun - b.fun, a.der - b.der)
    -(a::Dual, b::Real) = a-cdual(b)
    -(b::Real, a::Dual) = a-cdual(b)
    *(a::Dual, b::Dual) = Dual(a.fun * b.fun, a.fun*b.der+b.fun*a.der)
    *(b::Real, a::Dual) = a*cdual(b)
    *(a::Dual, b::Real) = a*cdual(b)
    /(a::Dual, b::Dual) = Dual(a.fun/b.fun, (a.der-(a.fun/b.fun)*b.der)/b.fun)
    /(a::Dual, b::Real) = a/cdual(b)
    /(b::Real, a::Dual) = cdual(b)/a
    ^(a::Dual, b::Real) = Dual(a.fun^b, b*a.fun^(b-1)*a.der)
    ==(a::Dual, b::Dual) = a.fun == b.fun && a.der == b.der
    
    function xdual(x0)
        Dual(x0,1)
    end
    
    cdual(x_0) = Dual(x_0,0)
end
