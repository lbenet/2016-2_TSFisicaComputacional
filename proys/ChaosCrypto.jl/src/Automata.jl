__precompile__(true)

"""
#Modulo AD (Lorenz.jl)

Diferenciación automática por medio de series de Taylor. 
Declaraciones del tipo, funciones que determinan cómo se realizan operaciones
manipulando series de Taylor. Base para trabajar con integración.
"""

module AD
    export Taylor, paso2,paso1
    export igualdad,logo,expo,seno,coseno

    type Taylor{T<:Number,S<:Number}
        ini :: T # Valor en el que se calcula la serie
        coef :: Array{S,1} # coeficiente normalizado de Taylor
    end
    Taylor(x0,a) = Taylor(promote(x0,a)...)

    function paso2(f,g)
    nocero = 1
    orden = length(f)
    h = zeros(f)
    total = 0
    while g[nocero] == 0
        nocero += 1
    end
    h[1] = f[nocero]/g[nocero]
    for k in nocero+1:orden
        for i in 1:k-nocero
            if length(h) < i || length(g) < k-i+1
                continue
            else
                total += h[i]*g[k-i+1]
            end
        end
        h[k-nocero+1] = (f[k] - total)/g[nocero]
    end
    h
    end

    function paso1(a,b)
    order=length(a)+length(b)
    c = zeros(order)
    total = 0
    i = 0
    for k in 1:order
        for i in 1:k
            if i >length(a) || k-i+1 > length(b)
                continue
            else
                total += a[i]*b[k-i+1]
            end
        end
        c[k] = total
        total = 0
    end
    c
    end

    function igualdad(a,b)
    orden = min([length(a),length(b)]...)
    i = 1
    while i <= orden
        if a[i] != b[i]
            return false
        end
        i+=1
    end
    return true
    end

    import Base: +, -, *, /, ==

    # Aqui se implementan los métodos necesarios para cada función

    +(a::Taylor, b::Taylor) = Taylor(a.ini,a.coef.+b.coef)
    -(a::Taylor, b::Taylor) = Taylor(a.ini,a.coef.-b.coef)
    *(a::Taylor, b::Taylor) = Taylor(a.ini,paso1(a.coef,b.coef))
    /(a::Taylor, b::Taylor) = Taylor(a.ini,paso2(a.coef,b.coef))
    ==(a::Taylor, b::Taylor) = igualdad(a.coef,b.coef)

    function logo(g,x0)
    L = zeros(g)
    L[1] = log(x0+1)
    total = 0
    for k in 0:length(g)-2
        total = 0
        for l in 0:k-1
            total += (l+1)*L[l+2]*g[k-l+1]
        end
        L[k+2] = ((k+1)*g[k+2]-total)/((k+1)*g[1])
    end
    Taylor(x0,L)
    end

    function expo(g,α,x0)
    P = zeros(g)
    P[1] = x0^(α)
    total = 0
    for k in 0:length(g)-2
        total = 0
        for l in 0:k-1
            total += α*((k-l+1)*P[l+1]*g[k-l+2])
            total -= (l+1)*P[l+2]*g[k-l+1]
        end
        total += α*P[k+1]*g[2]
        P[k+2] = total/((k+1)*g[1])
    end
    Taylor(x0,P)
    end

    function seno(g,x0)
    S = zeros(g)
    S[2] = cos(x0)
    total = 0
    totalidad = 0
    for k in 2:length(g)-1 # el término de inicio es truculento ya que si ponemos 1 sobreescribe S[2] = 1.
        total = 0
        totalidad = 0
        for l in 1:k-1
            for lp in 1:k-l
                total += lp*g[lp+1]*S[k-lp-l+1]
            end
            totalidad += l*g[l+1]*total/(k-l)
        end
        S[k+1] = -totalidad/k
    end
    Taylor(x0,S)
    end

    function coseno(g,x0)
    C = zeros(g)
    C[1] = cos(x0)
    total = 0
    totalidad = 0
    for k in 1:length(g)-1
        total = 0
        totalidad = 0
        for l in 1:k-1
            for lp in 1:k-l
                total += lp*g[lp+1]*C[k-lp-l+1]
            end
            totalidad += l*g[l+1]*total/(k-l)
        end
        C[k+1] = -totalidad/k
    end
    Taylor(x0,C)
    end

    import Base: log, cos, sin, ^

    cos(a::Taylor) = coseno(a.coef,a.ini)
    sin(a::Taylor) = seno(a.coef,a.ini)
    log(a::Taylor) = logo(a.coef,a.ini)
    ^(a::Taylor, α::Integer) = expo(a.coef,α,a.ini)
end
