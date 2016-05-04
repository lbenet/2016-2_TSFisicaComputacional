module SerieTaylor

#Importamos de base lo que vamos a necesitar
import Base: length, zero, one, convert, eltype, promote, promote_rule
import Base: log, sin, cos
import Base: +, -, *, /, ==, ^

export Taylor

#Declaramos el type Taylor
type Taylor{T<:Number}
    coeficientes :: Array{T,1}
    orden :: Int
    #Constructor interno
    function Taylor(coeficientes::Array{T,1}, orden::Int)  
        longitud_coeficientes= length(coeficientes)
        orden = max(orden, longitud_coeficientes-1)
        vector = zeros(T, orden+1)
        vector[1:longitud_coeficientes] = coeficientes[1:end]
        new(vector, orden)
    end
end

#Constructores externos
Taylor{T<:Number}(x::Taylor{T}, orden::Int) = Taylor{T}(x.coeficientes, orden)
Taylor{T<:Number}(x::Taylor{T}) = Taylor{T}(x.coeficientes, x.orden)
Taylor{T<:Number}(coeficientes::Array{T,1}, orden::Int) = Taylor{T}(coeficientes, orden)
Taylor{T<:Number}(coeficientes::Array{T,1}) = Taylor{T}(coeficientes, length(coeficientes)-1)
Taylor{T<:Number}(x::T, orden::Int) = Taylor{T}([x], orden)
Taylor{T<:Number}(x::T) = Taylor{T}([x], 0)

#Conversion y reglas de promocion
convert{T<:Number}(::Type{Taylor{T}}, a::Taylor) = Taylor(convert(Array{T,1}, a.coeficientes), a.orden)
convert{T<:Number, S<:Number}(::Type{Taylor{T}}, b::Array{S,1}) = Taylor(convert(Array{T,1},b))
convert{T<:Number, S<:Number}(::Type{Taylor{T}}, b::S) = Taylor([convert(T,b)], 0)
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{Taylor{S}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{Array{S,1}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Array{S,1}}, ::Type{Taylor{T}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{S}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{S}, ::Type{Taylor{T}}) = Taylor{promote_type(T, S)}

#Funcion auxiliar que nos declara el primer 'no cero'
function auxnocero{T<:Number}(a::Taylor{T})
    nocero::Int = a.orden+1
    for i in eachindex(a.coeficientes)
        if a.coeficientes[i] != zero(T)
            nocero = i-1
            break
        end
    end
    nocero
end

#Funcion auxiliar que nos arregla la forma de T y S
function auxforma{T<:Number, S<:Number}(a::Taylor{T}, b::Taylor{S})
    orden = max(a.orden, b.orden)
    a1, b1 = promote(a, b)
    return Taylor(a1, orden), Taylor(b1, orden), orden
end

#Tipo, longitud
eltype{T<:Number}(::Taylor{T}) = T
length{T<:Number}(a::Taylor{T}) = a.orden

#Cero y uno
zero{T<:Number}(a::Taylor{T}) = Taylor(zero(T), a.orden)
one{T<:Number}(a::Taylor{T}) = Taylor(one(T), a.orden)

#Igualdad
function ==(a::Taylor, b::Taylor)
    a1, b1, orden = auxforma(a, b)
    return a1.coeficientes == b1.coeficientes
end
==(a::Taylor, b::Number) = ==(a, Taylor(b, a.orden))
==(a::Number, b::Taylor) = ==(b, Taylor(a, b.orden))

#Suma y resta (con un poco de code generation)
for f in (:+, :-)
    @eval begin
        function ($f)(a::Taylor, b::Taylor)
            a1, b1, orden = auxforma(a, b)
            v = ($f)(a1.coeficientes, b1.coeficientes)
            return Taylor(v, orden)
        end
       ($f)(a::Taylor, b::Number) = ($f)(a, Taylor(b, a.orden))
       ($f)(a::Number, b::Taylor) = ($f)(Taylor(a, b.orden), b)
       ($f)(a::Taylor) = Taylor(($f)(a.coeficientes), a.orden)
    end
end

#Multiplicacion
function *(a::Taylor, b::Taylor)
    a1, b1, orden = auxforma(a, b)
    T = eltype(a1)
    coeficientes = zeros(T,orden+1)
    coeficientes[1] = a1.coeficientes[1] * b1.coeficientes[1]
    #Introducimos en el siguiente ciclo for a una 
    #funcion que declaramos justo despues,
    #el coeficiente homogeneo para la multiplicacion
    for k = 1:orden
        coeficientes[k+1] = mulHomogcoef(k, a1.coeficientes, b1.coeficientes)
    end
    Taylor(coeficientes, orden)
end

#Coeficiente homogeneo para multiplicacion
function mulHomogcoef{T<:Number}(kcoef::Int, ac::Array{T,1}, bc::Array{T,1})
    coefhomog = zero(T)
    for i = 0:kcoef
        coefhomog += ac[i+1] * bc[kcoef-i+1]
    end
    coefhomog
end
*(a::Taylor, b::Number) = Taylor(b*a.coeficientes, a.orden)
*(a::Number, b::Taylor) = Taylor(a*b.coeficientes, b.orden)

#Division
function /(a::Taylor, b::Taylor)
    a1, b1, orden = auxforma(a, b)
    ordLHopital, cLHopital = divlhopital(a1, b1) # orden y coeficiente L'Hopital
    T = typeof(cLHopital)
    v1 = convert(Array{T,1}, a1.coeficientes)
    v2 = convert(Array{T,1}, b1.coeficientes)
    coeficientes = zeros(T,orden+1)
    coeficientes[1] = cLHopital
    for k = ordLHopital+1:orden
        coeficientes[k-ordLHopital+1] = divHomogcoef(k, v1, v2, coeficientes, ordLHopital)
    end
    Taylor(coeficientes, orden)
end

#Funcion auxiliar que calcula orden L'Hopital, se asume a1 y b1 son del mismo orden
function divlhopital(a1::Taylor, b1::Taylor)
    a1nz = auxnocero(a1)
    b1nz = auxnocero(b1)
    ordLHopital = min(a1nz, b1nz)
    if ordLHopital > a1.orden
        ordLHopital = a1.orden
    end
    cLHopital = a1.coeficientes[ordLHopital+1] / b1.coeficientes[ordLHopital+1]
    auxsq = abs2(cLHopital)
    #Hacemos la prueba para ver si se puede aplicar L'Hopital
    if isinf(auxsq)
        info("Orden k=$(ordLHopital) => coeficientes[$(ordLHopital+1)]=$(cLHopital)")
        error("Division no define polinomio Taylor, o su primer coeficiente es infinito.\n")
    elseif isnan(auxsq)
        info("Orden k=$(ordLHopital) => coeficientes[$(ordLHopital+1)]=$(cLHopital)")
        error("No se puede aplicar L'Hopital.\n")
    elseif ordLHopital>0
        warn("Aplicando L'Hopital. Ultimos k=$(ordLHopital) coeficientes Taylor son 0.\n")
    end
    return ordLHopital, cLHopital
end

#Coeficiente Homogeneo para Division
function divHomogcoef{T<:Number}(kcoef::Integer, ac::Array{T,1}, bc::Array{T,1}, 
        coeficientes::Array{T,1}, ordLHopital::Integer)
    coefhomog = mulHomogcoef(kcoef, coeficientes, bc)
    coefhomog = (ac[kcoef+1]-coefhomog) / bc[ordLHopital+1]
    coefhomog
end
/(a::Taylor,b::Number) = Taylor(a.coeficientes/b, a.orden)
/(a::Number,b::Taylor) = Taylor([a], b.orden) / b

#Logaritmo
function log(a::Taylor)
    orden = a.orden
    l0nz = auxnocero(a)
    if auxnocero(a)>0
        error("No es posible expander log alrededor de 0.\n")
    end
    auxlog = log( a.coeficientes[1] )
    T = typeof(auxlog)
    ac = convert(Array{T,1}, a.coeficientes)
    coeficientes = zeros(T, orden+1)
    coeficientes[1] = auxlog
    for k = 1:orden
        coeficientes[k+1] = logHomogcoef(k, ac, coeficientes)
    end
    Taylor( coeficientes, orden )
end

#Coeficientes homogeneos para logaritmo
function logHomogcoef{T<:Number}(kcoef::Integer, ac::Array{T,1}, coeficientes::Array{T,1})
  coefhomog = zero(T)
  for i = 1:kcoef-1
    coefhomog += (kcoef-i) * ac[i+1] * coeficientes[kcoef-i+1]
  end
  coefhomog = (ac[kcoef+1] -coefhomog/kcoef) / ac[1]
  coefhomog
end

#Elevar al cuadrado
function square{T<:Number}(a::Taylor{T})
    orden = a.orden
    coeficientes = zeros(T,orden+1)
    coeficientes[1] = a.coeficientes[1]^2
    for k = 1:orden
        coeficientes[k+1] = squareHomogcoeff(k, a.coeficientes)
    end
    Taylor(coeficientes,orden)
end
#Coeficientes homogeneos para elevar al cuadrado
function squareHomogcoeff{T<:Number}(kcoef::Integer, ac::Array{T,1})
    coefhomog = zero(T)
    kodd = kcoef%2
    kend = div(kcoef - 2 + kodd, 2)
    for i = 0:kend
        coefhomog += ac[i+1]*ac[kcoef-i+1]
    end
    coefhomog = 2coefhomog
    if kodd == 0
        coefhomog += ac[div(kcoef,2)+1]^2
    end
    coefhomog
end

#Elevar a potencia entera cualquiera
function ^(a::Taylor, x::Integer)
    uno = one(a)
    if x < 0
        return uno / a^(-x)
    elseif x == 0
        return uno
    elseif x%2 == 0 # par
        if x == 2
            return square(a)
        else
            pow = div(x, 2)
            return square( a^pow )
        end
    else  # impar
        if x == 1
            return a
        else
            expon = div(x-1, 2)
            return a*square( a^expon )
        end
    end
end

#Sin y Cos
sin(a::Taylor) = sincos(a, "sin")
cos(a::Taylor) = sincos(a, "cos")
function sincos(a::Taylor, fun::AbstractString)
    orden = a.orden
    auxtrigo = sin( a.coeficientes[1] )
    T = typeof(auxtrigo)
    v = convert(Array{T,1}, a.coeficientes)
    sincoeficientes = zeros(T,orden+1)
    coscoeficientes = zeros(T,orden+1)
    sincoeficientes[1] = auxtrigo
    coscoeficientes[1] = cos( a.coeficientes[1] )
    for k = 1:orden
        sincoeficientes[k+1], coscoeficientes[k+1] = sincosHomogcoef(k, v, sincoeficientes, coscoeficientes)
    end
    if fun == "sin"
        return Taylor( sincoeficientes, orden )
    else
        return Taylor( coscoeficientes, orden )
    end
end

#Coeficientes homogeneos para cos y sin
function sincosHomogcoef{T<:Number}(kcoef::Integer, ac::Array{T,1}, 
        sincoeficientes::Array{T,1}, coscoeficientes::Array{T,1})
    sincoefhom = zero(T)
    coscoefhom = zero(T)
    for i = 1:kcoef
        number = i * ac[i+1]
        sincoefhom += number * coscoeficientes[kcoef-i+1]
        coscoefhom -= number * sincoeficientes[kcoef-i+1]
    end
    sincoefhom = sincoefhom/kcoef
    coscoefhom = coscoefhom/kcoef
    return sincoefhom, coscoefhom
end
    




