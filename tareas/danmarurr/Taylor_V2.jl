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


#Importamos todas las funciones para las cuales queremos definir su operacion 
import Base: ^, exp, sqrt, cbrt, sin, cos, tan, cot, sec, csc, sinh, cosh, tanh, coth, sech, csch,
asin,   acos,   atan,   acot,   asec,   acsc,
asinh,  acosh,  atanh,  acoth,  asech,  acsch

#Vector con todos los símbolos asociados a las funciones trigométricas, exponenciales, etc. y sus respectivas
#derivadas
Vec_Func = [(:sin, :cos), (:cos, :(x -> -sin(x))), (:tan, :(x -> (sec(x))^2)), (:cot, :(x -> -(csc(x))^2)), 
    (:sec, :(x -> sec(x)*tan(x))), (:csc, :(x -> -csc(x)*cot(x))), (:sinh, :cosh), (:cosh, :sinh), 
    (:tanh, :(x -> (sech(x))^2)), (:coth, :(x -> -(csch(x))^2)), (:asin, :(x -> 1/sqrt(1-x^2))), 
    (:acos, :(x -> -1/sqrt(1-x^2))), (:atan, :(x -> 1/(1+x^2))), (:acot, :(x -> -1/(1+x^2))),
    (:asec, :(x -> 1/(sqrt(x^2-1)*x))), (:acsc, :(x -> -1/(sqrt(x^2-1)*x))), (:asinh, :(x -> 1/sqrt(1+x^2))),
    (:acosh, :(x -> 1/sqrt(x^2-1))), (:atanh, :(x -> 1/(1-x^2))), (:acoth, :(x -> 1/(1-x^2))),
    (:asech, :(x -> 1/(x*sqrt(1-x^2)))), (:acsch, :(x -> -1/(x*sqrt(1+x^2)))), (:sqrt, :(x -> 1/(2*sqrt(x)))), (:exp, :exp),
    (:cbrt, :(x -> 1/(3*x^(2/3))))]

#Casos especiales: los logarítmos y a^x
log(a::Dual) = Dual(log(a.fun), a.der/a.fun)
log{T<:Real}(b::T, a::Dual) = Dual(log(b,a.fun), a.der/(log(b)*a.fun))
^{T<:Real}(b::T, a::Dual) = Dual(b^a.fun, a.der*log(b)*b^a.fun)

#Loop para crear los nuevos métodos a partir de los símbolos
for r in 1:length(Vec_Func)
    fn = Vec_Func[r][1] #El primer símbolo está asociado a la función
    der = Vec_Func[r][2] #El segundo símbolo está asociado a la derivada de la función
    ex = quote #Creamos una nueva expresión
        function ($fn)(a::Dual) #Definimos fn(a::Dual)
            fun = ($fn)(a.fun)
            derv = ($der)(a.fun)
            return Dual(fun, derv*a.der) #Aplicamos la regla de la cadena
        end
    end
    @eval $ex #Evaluamos la expresión para crear el método
end

# ***********************************************************
# DERVIADAS DE ORDEN SUPERIOR
# ***********************************************************

"""Definición de polinomios de Taylor, donde
...
"""
type Taylor{T<:Number}
    coef :: Array{T,1}     #Arreglo de coeficientes
    order :: Int64         # Orden del polinimio
end

#Constructores para el tipo Taylor
Taylor{T<:Number}(v::Array{T,1}) = Taylor(v, length(v))   #Recibe un arreglo de coeficientes
Taylor(N::Int64) = Taylor(zeros(N), N) #Recibe solo el orden y llena con ceros.


function upgrade(a::Taylor, b::Taylor)
    order_a = a.order
    order_b = b.order
    dif = order_a - order_b    #Calcula la diferencia entre los ordenes de dos Taylor
    
    if dif >= 0     
        for i in 1:abs(dif)
            push!(b.coef, 0)  #Si el orden de b es menor, agrega tantos ceros como hagan falta para igualar a a.
        end
        b.order = a.order
    else
        for i in 1:abs(dif)
            push!(a.coef, 0)  # Lo mismo si a es de orden menor.
        end
        a.order = b.order
    end
    return a,b
        
    
    
end

# Para la suma y la resta usamos metaprogramming

for fn = (:+, :-)
    ex = quote
        function ($fn)(a::Taylor, b::Taylor)
            c,d = upgrade(a,b)     # Hace un upgrade, asi no es tragico si los Taylor tienen orden diferente.
            s = Taylor(c.order)
            for i in 1:s.order
                s.coef[i] = ($fn)(c.coef[i], d.coef[i])
            end
            return s
        end
    end
        @eval $ex
end

function /(a::Taylor, b::Taylor)
    f, g = upgrade(a,b)
    k = f.order
    l = g.order
    r = 0
    for j in 1:l
        if g.coef[j] != 0 
            r = j
            break
        end
    end
    @assert r != 0 "Requerimos que g sea distinto de 0"
    h = Taylor(k)
    #@show g
    h.coef[1] = f.coef[r] / g.coef[r]
    
    for j in (r + 1):k
        #@show j
        f_j = f.coef[j]
        suma = 0
        for i in 1:(j - r)
            #@show i
            h_i = h.coef[i]
            g_i = g.coef[j - i + 1]
            suma += h_i * g_i
        end
        h.coef[j + 1 - r] = (f_j - suma)/g.coef[r]
    end
    return h
end

function *(a::Taylor, b::Taylor)
    f,g = upgrade(a,b)
    h = Taylor(f.order)
    
    for k in 1:f.order
        s = 0
        for i in 1:k
            
            h.coef[k] += f.coef[i]*g.coef[(k+1) - i]
            #println(k, i, "  ", s)
            
        
        end
        
    end
    return h
            
end

function *{T<:Number}(f::Taylor, α::T)
    h = Taylor(zeros(typeof(promote(α, f.coef[1])[1]), f.order))
    
    for k in 1:f.order
        h.coef[k] = α * f.coef[k]        
    end
    return h
            
end
*{T<:Number}(α::T, f::Taylor) = *(f, α)

function ==(a::Taylor, b::Taylor)
    f,g = upgrade(a,b)
    r = true
    for k in 1:f.order
        if g.coef[k] != f.coef[k]
            r = false
            break
        end
    end
    return r
end

# *************** FUNCIONES ESPECIALES **********************


function log(g::Taylor)
    L = Taylor(g.order)
    r = 0
    
    for i in 1:g.order
        if g.coef[i] != 0
            r = i
            
            break
        end
    end
    @assert r != 0 "Necesitamos un polinomio distinto de cero"
    #@show r
    g_r = g.coef[r]
    L.coef[r] = log(g_r)
    
    for k in (r + 1):g.order
        
        k != g.order ? g_k = g.coef[k + 1] : g_k = 0
        suma = 0
        for j in(r+1):k
            g_j = g.coef[j]
            L_j = L.coef[k - j + 1]
            suma += (k - j)*g_j*L_j
        end
        L.coef[k - r + 1] = (k*g_k - suma)/(g_r*(k - r))
    end
    return L
end

function exp(a::Taylor)
    E = Taylor(zeros(typeof(exp(a.coef[1])),a.order))     #Esta exponencial permite coeficientes complejos
    r = 0
    for i in 1:a.order #Verificamos que el Taylor proporcionado sea distinto de cero
        if a.coef[i] != 0
            r = i
            break
        end
    end
    @assert r != 0 "Necesitamos un polinomio distinto de cero"
    E.coef[1] = exp(a.coef[1])
    for k in 2:a.order
        suma = 0
        for j in 1:k-1
            suma += E.coef[j] * (k - j) * a.coef[k-j+1]
        end
        E.coef[k] = suma / (k-1)
    end
    E
end

import Base: real, imag

function real(f::Taylor)
    h = Taylor(zeros(typeof(real(f.coef[1])),f.order))
    
    for k in 1:f.order
        h.coef[k] = real(f.coef[k])
    end
    return h
end

function imag(f::Taylor)
    h = Taylor(zeros(typeof(real(f.coef[1])),f.order))
    
    for k in 1:f.order
        h.coef[k] = imag(f.coef[k])
    end
    return h
end

function cos(f::Taylor)
    t1 = im*f
    ex = exp(t1)
    return real(ex)
end

function sin(f::Taylor)
    t1 = im*f
    ex = exp(t1)
    return imag(ex)
end

^{T<:Number}(a::Taylor, n::T) = exp(n*log(a))


end

