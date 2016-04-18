type Taylor{T<:Number}
    
    coeficientes :: Array{T,1}
    orden :: T
    
    #Constructor interno
    function Taylor(coeficientes::Array{T,1}, orden::T)
        
        longitud_coeficientes= length(coeficientes)
        orden = max(orden, longitud_coeficientes-1)
        
        vector = zeros(T, orden+1)
        
        vector[1:longitud_coeficientes] = coeficientes[1:end]
        new(vector, orden)
        
    end
end

#Constructor externo
Taylor{T<:Number}(x::Taylor{T}, orden::Integer) = Taylor{T}(x.coeficientes, orden)
#Taylor{T<:Number}(x::Taylor{T}) = Taylor{T}(x.coeficientes, x.orden)
Taylor{T<:Number}(x::Taylor{T}) = x
Taylor{T<:Number}(coeficientes::Array{T,1}, orden::Integer) = Taylor{T}(coeficientes, orden)
Taylor{T<:Number}(coeficientes::Array{T,1}) = Taylor{T}(coeficientes, length(coeficientes)-1)
Taylor{T<:Number}(x::T, orden::Integer) = Taylor{T}([x], orden)
Taylor{T<:Number}(x::T) = Taylor{T}([x], 0)

#Tipo, longitud
eltype{T<:Number}(::Taylor{T}) = T
length{T<:Number}(a::Taylor{T}) = a.orden

#Conversion y reglas promocion
convert{T<:Number}(::Type{Taylor{T}}, a::Taylor) = Taylor(convert(Array{T,1}, a.coeficientes), a.orden)
convert{T<:Number, S<:Number}(::Type{Taylor{T}}, b::Array{S,1}) = Taylor(convert(Array{T,1},b))
convert{T<:Number, S<:Number}(::Type{Taylor{T}}, b::S) = Taylor([convert(T,b)], 0)

promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{Taylor{S}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{Array{S,1}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Array{S,1}}, ::Type{Taylor{T}}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{Taylor{T}}, ::Type{S}) = Taylor{promote_type(T, S)}
promote_rule{T<:Number, S<:Number}(::Type{S}, ::Type{Taylor{T}}) = Taylor{promote_type(T, S)}

# #Funcion auxiliar
function firstnonzero{T<:Number}(a::Taylor{T})
    orden = a.orden
    nonzero = orden+1
    z = zero(T)
    for i = 1:orden+1
        if a.coeficientes[i] != z
            nonzero = i-1
            break
        end
    end
    nonzero
end
function fixshape{T<:Number, S<:Number}(a::Taylor{T}, b::Taylor{S})
    orden = max(a.orden, b.orden)
    a1, b1 = promote(a, b)
    return Taylor(a1, orden), Taylor(b1, orden), orden
end

#Zero y uno
zero{T<:Number}(a::Taylor{T}) = Taylor(zero(T), a.orden)
one{T<:Number}(a::Taylor{T}) = Taylor(one(T), a.orden)

import Base: +, -, *, /, ==

#Igualdad
function ==(a::Taylor, b::Taylor)
    a1, b1, orden = fixshape(a, b)
    return a1.coeficientes == b1.coeficientes
end
==(a::Taylor, b::Number) = ==(a, Taylor(b, a.orden))
==(a::Number, b::Taylor) = ==(b, Taylor(a, b.orden))

#Suma y resta
for f in (:+, :-)
    @eval begin
        function ($f)(a::Taylor, b::Taylor)
            a1, b1, orden = fixshape(a, b)
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
    a1, b1, orden = fixshape(a, b)
    T = eltype(a1)
    coeficientes = zeros(T,orden+1)
    coeficientes[1] = a1.coeficientes[1] * b1.coeficientes[1]
    for k = 1:orden
        coeficientes[k+1] = mulHomogcoef(k, a1.coeficientes, b1.coeficientes)
    end
    Taylor(coeficientes, orden)
end

#Coeficiente homogeneo para multiplicacion
function mulHomogcoef{T<:Number}(kcoef::Integer, ac::Array{T,1}, bc::Array{T,1})
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
    a1, b1, orden = fixshape(a, b)
    ordLHopital, cLHopital = divlhopital(a1, b1) # L'Hôpital orden y coeficiente
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

# Coeficiente Homogeneo para Division
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
    l0nz = firstnonzero(a)
    if firstnonzero(a)>0
        error("No es posible expander log alrededor de 0.\n")
    end
    aux = log( a.coeficientes[1] )
    T = typeof(aux)
    ac = convert(Array{T,1}, a.coeficientes)
    coeficientes = zeros(T, orden+1)
    coeficientes[1] = aux
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

#Exponente
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
function sincos(a::Taylor, fun::String)
    orden = a.orden
    aux = sin( a.coeficientes[1] )
    T = typeof(aux)
    v = convert(Array{T,1}, a.coeficientes)
    sincoeficientes = zeros(T,orden+1)
    coscoeficientes = zeros(T,orden+1)
    sincoeficientes[1] = aux
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








