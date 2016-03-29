"""
Modulo de Diferenciación Automática

https://en.wikipedia.org/wiki/Automatic_differentiation

Duales(Real,Real)
xdual(Real)
Base: +, -, *, /, ^

Autor: Carlos López
22/03/2016
"""
__precompile__()

module AD
    import Base: +, -, *, /, ^

    type Dual
        fun::Real
        der::Real
    end
    Dual(a::Real) = Dual(a,0)
    Dual(a::Irrational,b::Irrational) = Dual(convert(Float64,a),convert(Float64,b))

    xdual(x0)=Dual(x0,1)
    

    import Base: +, -, *, /, ^, == 
    +(x::Dual,y::Dual) =Dual( x.fun+y.fun , x.der+y.der )
    -(x::Dual,y::Dual) =Dual( x.fun-y.fun , x.der-y.der )
    *(x::Dual,y::Dual) =Dual( x.fun*y.fun , (x.fun*y.der)+(y.fun*x.der) )
    /(x::Dual,y::Dual) =Dual( x.fun/y.fun , (x.der-((x.fun/y.fun)*y.der))/y.fun )

    +(x::Dual,y::Real) =Dual( x.fun+y , x.der )
    -(x::Dual,y::Real) =Dual( x.fun-y , x.der )
    *(x::Dual,y::Real) =x*convert(Dual,y)
    /(x::Dual,y::Real) =x/convert(Dual,y)

    +(y::Real,x::Dual) =Dual( x.fun+y , x.der )
    -(y::Real,x::Dual) =Dual( x.fun-y , x.der )
    *(y::Real,x::Dual) =x*convert(Dual,y)
    /(y::Real,x::Dual) =x/convert(Dual,y)

    ^(x::Dual,a::Integer) =Dual( x.fun^a , a*(x.fun^(a-1))*x.der )
    ^(x::Dual,a::Real) =Dual( x.fun^a , a*(x.fun^(a-1))*x.der )

    +(x::Dual) =Dual( +x.fun , +x.der )
    -(x::Dual) =Dual( -x.fun , -x.der )

    ==(x::Dual,y::Dual)=(x.fun,x.der)==(y.fun,y.der)


    import Base: convert, promote_rule
    convert(::Type{Dual},x::Real) = Dual(x,0)
    promote_rule(::Type{Dual}, ::Type{Int} ) = Dual
    promote_rule(::Type{Dual}, ::Type{Float64} ) = Dual


    import Base: print,show
    show(io::IO,x::Dual)=print(io, x.fun," + ",x.der,"e")
    print(io::IO,x::Dual)=show(io,x)

    export Dual, xdual

end