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
    import Base: +, -, *, /, ^, ==

    type Dual
        fun::Real
        der::Real
    end
    Dual(a::Real) = Dual(a,0)

    Dual(a::Irrational,b::Irrational) = Dual(convert(Float64,a),convert(Float64,b))
    Dual(a::Irrational,b::Real) = Dual(convert(Float64,a),convert(Float64,b))
    Dual(a::Real,b::Irrational) = Dual(convert(Float64,a),convert(Float64,b))

    xdual(x0)=Dual(x0,1)
    
    
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
    /(y::Real,x::Dual) =convert(Dual,y)/x

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
    show(io::IO,x::Dual)=print(io, x.fun," + ",x.der,"ε")
    print(io::IO,x::Dual)=show(io,x)

    export Dual, xdual


    derivadas = (
    [:exp,:(exp(x.fun))],
    [:sqrt,:(1/(2*sqrt(x.fun)))],
    [:log,:(1/x.fun)],
    [:sin,:(cos(x.fun))],
    [:cos,:(-sin(x.fun))],
    [:tan,:(sec(x.fun)^2)],
    [:cot,:(-csc(x.fun)^2)],
    [:sec,:(sec(x.fun)*tan(x.fun))],
    [:csc,:(-csc(x.fun)*cot(x.fun))],
    [:sinh,:(cosh(x.fun))],
    [:cosh,:(sinh(x.fun))],
    [:tanh,:(sec(x.fun)^2)],
    [:coth,:(-csch(x.fun)^2)],
    [:sech,:(sech(x.fun)*tanh(x.fun))],
    [:csch,:(-csch(x.fun)*coth(x.fun))],
    [:asin,:(1/sqrt(1-(x.fun^2)))],
    [:acos,:(-1/sqrt(1-(x.fun^2)))],
    [:atan,:(1/sqrt(1+(x.fun^2)))],
    [:acot,:(-1/(1+(x.fun^2)))],
    [:asec,:(1/(abs(x.fun)*sqrt((x.fun^2)-1)))],
    [:acsc,:(-1/(abs(x.fun)*sqrt((x.fun^2)-1)))],
    [:asinh,:(1/sqrt((x.fun^2)+1))],
    [:acosh,:(1/sqrt((x.fun^2)-1))],
    [:atanh,:(1/sqrt(1-(x.fun^2)))],
    [:acoth,:(1/sqrt(1-(x.fun^2)))],
    [:asech,:(1/(1-(x.fun^2)))],
    [:acsch,:(-1/(abs(x.fun)*sqrt((x.fun^2)+1)))],
    [:sinc,:(cosc(x.fun))],
    )

    for op in derivadas
      eval(quote
            import Base:$(op[1])
            ($(op[1]))(x::Dual) = Dual( $(op[1])(x.fun) , x.der*$(op[2]) ) 
      end)
    end

end