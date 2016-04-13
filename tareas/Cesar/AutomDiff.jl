#22/03/16

#Este módulo implementa la definición del tipo Taylor, así como la de sus operaciones.

#Rodríguez Rosenblueth César Daniel
#Santiago Santos Eva Yazmín
#Hernández de la Vega Alejandro

__precompile__(true)

module AD
    import Base: +, -, *, /, ^, ==

export Taylor

#Se define el tipo Taylor:
type Taylor{D <: Number}
    coef :: Array{D,1}
    order :: Int
    

  #Definimos una función para que considere el coeficiente y el orden de la serie de Taylor
    function Taylor(coef :: Array{D,1}, order :: Int)
        l = length(coef) #Definimos la longitud del Taylor
        order = max(order,l-1) 
        order == l-1 && return new(coef, order)
        resize!(coef, order+1)
        for i = l+1:order+1
            coef[i]=zero(D)
        end
        new(coef, order)
    end
    
end

Taylor{D<:Number}(x::Taylor{D}, order::Int) = Taylor{D}(x.coef, order) # Un taylor siempre será un taylor
Taylor{D<:Number}(x::Taylor{D}) = x #El taylor de un taylor, es un taylor
Taylor{D<:Number}(coef::Array{D,1}, order::Int) = Taylor{D}(coef, order) # La entrada más sencilla
Taylor{D<:Number}(coef::Array{D,1}) = Taylor{D}(coef, length(coef)-1) # El taylor de un arreglo será un taylor
Taylor{D<:Number}(x::D, order::Int) = Taylor{D}([x], order) # De una variable
Taylor{D<:Number}(x::D) = Taylor{D}([x], 0) # De una constante





import Base: +, -, *, /, ^, ==

#Definimos con el mínimo para no tener que agregar ceros. De por sí por como definimos Taylor manualmente puedes especificar 
# el orden necesario.
function +(a::Taylor,b::Taylor)
    m = minimum([a.order,b.order])+1
    c = zeros(m)
    for k in range(1,m)
        c[k]=a.coef[k]+b.coef[k]
    end
    return Taylor(c)
end

function -(a::Taylor,b::Taylor)
    m = minimum([a.order,b.order])+1
    c = zeros(m)
    for k in range(1,m)
        c[k]=a.coef[k]-b.coef[k]
    end
    return Taylor(c)
end

function *(a::Taylor,b::Taylor)
    m = minimum([a.order,b.order])+1
    c = zeros(m)
    for k in range(1,m)
        for i in range(1,k)
            c[k]+=a.coef[i]*b.coef[k-i+1]
        end
    end
    return Taylor(c)
end

function /(a::Taylor,b::Taylor)
    m = minimum([a.order,b.order])+1
    c = zeros(m)
    for k in range(1,m)
        s = 0
        for i in range(1,k-1)
            s += c[i]*b.coef[k-i+1]
        end
        c[k] = (a.coef[k]-s)/b.coef[1]
    end
    return Taylor(c)
end

function /(a::Taylor,b::Taylor)
    for k in range(1,b.order+1)
        global l
        if b.coef[k] != 0
            l=k
            break
        end
    end
    m = maximum([a.order,l])
    c = zeros(m)
    a = Taylor(a.coef,m)
    b = Taylor(b.coef,m)
    for k in range(1,m-l+1)
        s = 0
        for i in range(1,k-l)
            s += c[i]*b.coef[k+l-i]
        end
        c[k] = (a.coef[k+l-1]-s)/b.coef[l]
    end
    return Taylor(c)
end

function ==(a::Taylor,b::Taylor)
    m = maximum([a.order,b.order])
    #Agregamos ceros
    a = Taylor(a.coef,m)
    b = Taylor(b.coef,m)
    for k in range(1,m+1)
        if a.coef[k] != b.coef[k]
            return false
        end
    end
    return true
end

end
