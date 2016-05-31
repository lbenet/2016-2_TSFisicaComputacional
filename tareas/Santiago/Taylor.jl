__precompile__(true)

module TS
	import Base: +, -, *, /, ^, ==
	import Base: exp, log, sin, cos
	import Base.promote

	export Taylor
"""
Definici�n de polinomios de Taylor, consta de dos partes: 'taylor_vec' es un Array{T<:Number} que contiene los
coeficientes ordenados del polinomio y 'gr' es el grado del polinomio de tal forma que
Taylor{Number}([a?,a?,a?,..,a_n],m) representa el polinomio a? + a?x + a?x� + ��� + a?x? + ... + 0x?
...
"""
type Taylor{T<:Number}
    # Aqu� se declara taylor_vec que es un Array unidimensional de tipo T.
    taylor_vec::Array{T, 1}
    # gr es un entero que corresponde al grado m�ximo de la serie.
    gr :: Int

    # Aqu� hay un constructor interno que lo que va a hacer es tomar el
    # arreglo taylor_vec y ver si tiene un tama�o mayor o menor que gr
    # que es donde se trunca la serie.
    function Taylor(taylor_vec::Array{T, 1},gr::Int)
        if length(taylor_vec) < gr + 1
            # Si hay menos de gr entradas, se utiliza vcat para meter ceros
            # restantes en taylor_vec y luego se reconstruye.
            taylor_vec = vcat(taylor_vec, zeros(T, gr + 1 - length(taylor_vec)))
            new(taylor_vec,gr)
        else
            # En caso contrario se seleccionan las primeras gr + 1 entradas.
            taylor_vec = taylor_vec[1:gr+1]
            new(taylor_vec,gr)
        end
    end
end

# Para no tener que especificar el tipo de los coeficientes:
Taylor{T<:Number}(taylor_vec::Array{T, 1},gr::Int) = Taylor{T}(taylor_vec,gr)

# Si no queremos especificar el grado de los Taylor:
Taylor{T<:Number}(a::Array{T,1}) = Taylor(a, length(a)-1)

# Adicionalmente se declara un m�todo para tener Taylor para un escalar.
Taylor{T<:Number}(a::T) = Taylor([a])
Taylor{T<:Number}(a::T, gr::Int) = Taylor([a],gr)

# M�todo para la promoci�n de dos Taylors y as� poder actuar
# estructuras aritm�ticas:
function promote(a::Taylor, b::Taylor)
    gr = maximum([a.gr, b.gr])
    return Taylor(a.taylor_vec,gr), Taylor(b.taylor_vec, gr)
end

# La suma entre los dos Taylors y entre un escalar.
function +(a::Taylor, b::Taylor)
    a, b = promote(a,b)
    return Taylor(a.taylor_vec + b.taylor_vec, a.gr)
end
+(a::Taylor, b::Number) = a + Taylor(b)
+(a::Number, b::Taylor) = Taylor(a) + b

# La resta se implementa igualmente.
function -(a::Taylor, b::Taylor)
    a, b = promote(a,b)
    return Taylor(a.taylor_vec - b.taylor_vec, a.gr)
end
-(a::Taylor, b::Number) = a - Taylor(b)
-(a::Number, b::Taylor) = Taylor(a) - b

# Los operadores unitarios y la comparaci�n.
+(a::Taylor) = a
-(a::Taylor) = Taylor(-a.taylor_vec,a.gr)
function ==(a::Taylor, b::Taylor)
    a, b = promote(a,b)
    return a.taylor_vec == b.taylor_vec
end

# La mutiplicaci�n se implementa usando dos for's. Primero se extienden
# Los arreglos porque la mutiplicaci�n de un polinomio de grado a.gr por
# uno de grado b.gr es un polinomio de grado a.gr + b.gr
# Una vez hecha la extensi�n, el primero for que va
# desde 1, length(a.taylor_vec) construye el arreglo entero de
# elementos. El segundo for anidado es la suma que se da en la teor�a
# de arriba. S�lo hay que tener cuidado porque los �ndices empiezan en
# uno y no en cero como en la teor�a.
function *(a::Taylor, b::Taylor)
    a, b = Taylor(a.taylor_vec,a.gr + b.gr), Taylor(b.taylor_vec, a.gr + b.gr)
    return Taylor([sum([a.taylor_vec[i]*b.taylor_vec[k-i+1]
            for i in range(1,k)]) for k in range(1,a.gr+1)],a.gr)
end
*(a::Taylor, b::Number) = a*Taylor(b)
*(a::Number, b::Taylor) = Taylor(a)*b

# La divisi�n resulta ser un poco m�s complicada:
function /(a::Taylor, b::Taylor)
    # Lo primero que se hace es ver el orden de la divisi�n:
    m = a.gr - b.gr
    # Si m es negativo se avienta un error.
    if m < 0
        return error("Se est� dividiendo un polinomio de orden k entre otro con orden mayor a k.")
    end
    # Ahora vamos a promover los Taylor a un mismo tama�o:
    a, b = promote(a,b)
    # Se inicializa un arreglo lleno de ceros enteros. Este arreglo despu�s se
    # va a cortar al tama�o m +1.
    div_vec = zeros(Int, a.gr+1)
    # Si el primer elemento de b es cero, nuestro algoritmo, tal cual como
    # est� escrito no funciona. Debe empezar con el primer par a.taylor_vec[i]
    # b.taylor_vec[i] con b.taylor_vec[i] distinto de cero. Esto se hace con
    # un loop sobre los elementos de b. j va a ser un contador que nos diga
    # d�nde es el coeficiente distinto de cero.
    j = 1
    for i in range(1,length(b.taylor_vec))
        if b.taylor_vec[i] != 0
            div_vec[1] = a.taylor_vec[i]/b.taylor_vec[i]
            j = i
            break
        # Si llegamos al final del arreglo sin encontrar un coeficiente de b
        # distinto de cero, regresamos un error.
        elseif (i==length(b.taylor_vec)) & (b.taylor_vec[i] == 0)
            return error("El polinomio divisor es igual a cero")
        end
    end
    # Habiendo obtenido la j que es donde el divisor empieza a ser distinto de cero,
    # debemos comprobar que todos los coeficientes de a anteriores a j sean cero. De
    # otra manera tendr�amos una divisi�n de un t�rmino x^k/x^j con k<j, lo cual no
    # tiene una expansi�n en Taylor.
    for i in range(1,j-1)
        if a.taylor_vec[i] != 0
            return error("Polinomio con polo en la expansi�n de Taylor.")
        end
    end
    # Ahora podemos llenar los elementos que falta. Nuestro rango empieza
    # en 2 pero no necesariamente termina en el �ndice correspondiente al tama�o
    # por la parte anterior que vimos que el primer elemento no necesariamente
    # se calcul� con �ndices 1,1. Para eso usamos el contador j
    for k in range(j+1, length(div_vec)-j)
        # sd ser� la suma que aparece en la ecuaci�n, inicializada como cero.
        sd = 0
        # Con el siguiente for se ejecuta la suma utilizando los coeficientes
        # de div_vec que ya se generaron.
        for i in range(1,k-j)
            sd = sd + (div_vec[i])*b.taylor_vec[k-i+1]
        end
        div_vec[1+k-j] = (1/b.taylor_vec[j])*(a.taylor_vec[k]-sd)
    end
    # Se implementa el orden correcto m
    return Taylor(div_vec,m)
end
/(a::Taylor, b::Number) = a/Taylor(b)
/(a::Number, b::Taylor) = Taylor(a)/b

# Esta es la funci�n exponencial aplicada a estructuras Taylor.
function exp(a::Taylor,gr::Int)
    # Este primer if est� armado para que no se regrese un error de
    # inexactitud de meter un entero a una exponencial que trabaja
    # m�nimo con n�meros flotantes.
    if eltype(a.taylor_vec) <: Integer
        a = Taylor(convert(Array{Float64,1}, a.taylor_vec))
    end
    # Inicializamos el Taylor e de grado gr.
    e = Taylor(0.,gr)
    # Promovemos los arreglos a y e.
    a, e = promote(a,e)
    # Ahora la condici�n inicial sobre el arreglo.
    e.taylor_vec[1] = exp(a.taylor_vec[1])
    # El siguiente for es la relaci�n de recurrencia hallada que
    # usa los valores E_k previamente obtenidos.
    for k in range(2,gr)
        e.taylor_vec[k] = (1/(k-1))*sum([(k-j)*e.taylor_vec[j]*a.taylor_vec[k-j+1] for j in range(1,k)])
    end
    # Se regresa la estructura de Taylor con coeficientes calculados.
    return e
end

# Implementamos la funci�n sin poner el orden expl�cito.
exp(a::Taylor) = exp(a,a.gr)


# Esta es la funci�n logaritmo aplicada a estructuras Taylor.
function log(a::Taylor,gr::Int)
    # Este primer if est� armado para que no se regrese un error de
    # inexactitud de meter un entero a un logaritmo que trabaja
    # m�nimo con n�meros flotantes.
    if eltype(a.taylor_vec) <: Integer
        a = Taylor(convert(Array{Float64,1}, a.taylor_vec))
    end
    # Debemos revisar que no se calcule el logaritmo de cero.
    if a.taylor_vec[1] == 0
        return error("No se puede calcular el logaritmo de cero.")
    end
    # Es posible calcular el logaritmo de n�meros negativos pero tendr�amos
    # que convertir el arreglo a n�meros complejos:
    if (eltype(a.taylor_vec[1]) <: Real) & (a.taylor_vec[1] < 0)
        a = Taylor(convert(Array{Complex{Float64},1}, a.taylor_vec))
    end
    # Inicializamos el Taylor l con la condici�n inicial:
    l = Taylor(log(a.taylor_vec[1]),gr)
    if a.gr == 0
        return l
    end
    l.taylor_vec[2] = a.taylor_vec[2]/a.taylor_vec[1]
    a, l = promote(a,l)
    # El siguiente for es la relaci�n de recurrencia hallada que
    # es b�sicamente el de la divisi�n de dos series de Taylor.
    for k in range(2,gr-1)
        l.taylor_vec[k+1] = (1/a.taylor_vec[1])*(a.taylor_vec[k+1]-(1/k)*sum([i*l.taylor_vec[i+1]*a.taylor_vec[k-i+1] for i in range(1,k-1)]))
    end
    # Se regresa la estructura de Taylor con coeficientes calculados.
    return l
end

log(a::Taylor) = log(a,a.gr)


# Esta es la funci�n potencia aplicada a Taylors
# Primero Julia nos pide que definamos para un exponente entero
# Por defautlt, si ex es entero se regresa un Taylor de orden ex veces
# el Taylor original.
function ^(a::Taylor,ex::Integer)
    if ex < 0
        return error("No se puede exponenciar la serie negativamente")
    end
    if a.gr == 0
        return Taylor(a.taylor_vec[1]^ex)
    end
    # Ahora encontramos donde empieza el arreglo que se está exponenciando,
    # análogo a la división.
    j = 1
    for k in range(1,a.gr+1)
        if a.taylor_vec[k] != 0
            j = k
            break
        elseif k==a.gr+1
            # Si el arreglo es un cero se regresa un cero.
            return Taylor(0*a.taylor_vec[1])
        end
    end
    # Factorizamos el arreglo:
    a = Taylor(a.taylor_vec[j:a.gr+1])
    p = Taylor(a.taylor_vec[1]^ex,a.gr*ex)
    a,p = promote(a,p)
    # El siguiente for es la relación de recurrencia hallada que
    # es básicamente el de la división de dos series de Taylor.
    for k in range(1,p.gr)
        p.taylor_vec[k+1] = (ex/k)*(1/a.taylor_vec[1])*(sum([p.taylor_vec[i]*(k+1-i)*a.taylor_vec[k-i+2] for i in range(1,k)])
        -sum([(i/ex)*p.taylor_vec[i+1]*a.taylor_vec[k-i+1] for i in range(1,k-1)]))
    end
    return Taylor(vcat(zeros(eltype(p.taylor_vec),(j-1)*ex),p.taylor_vec))
end

function ^(a::Taylor,ex::Number)
    if a.gr == 0
        return Taylor(a.taylor_vec[1]^ex)
    end
    # Ahora encontramos donde empieza el arreglo que se está exponenciando,
    # análogo a la división.
    if a.taylor_vec[1] ==0
        return error("No se puede exponenciar a una pontenica no-entera si el término constante es cero.")
    end
    p = Taylor(a.taylor_vec[1]^ex,a.gr)
    # El siguiente for es la relación de recurrencia hallada que
    # es básicamente el de la división de dos series de Taylor.
    for k in range(1,p.gr)
        p.taylor_vec[k+1] = (ex/k)*(1/a.taylor_vec[1])*(sum([p.taylor_vec[i]*(k+1-i)*a.taylor_vec[k-i+2] for i in range(1,k)])
        -sum([(i/ex)*p.taylor_vec[i+1]*a.taylor_vec[k-i+1] for i in range(1,k-1)]))
    end
    return p
end

^(a::Taylor,b::Number,c::Int) = ^(Taylor(a.taylor_vec,c),b)

function sin(a::Taylor,gr::Int)
    # Este primer if est� armado para que no se regrese un error de
    # inexactitud de meter un entero a una exponencial que trabaja
    # m�nimo con n�meros flotantes.
    if eltype(a.taylor_vec) <: Integer
        a = Taylor(convert(Array{Float64,1}, a.taylor_vec))
    end
    # Inicializamos los Taylor s,c de grado gr.
    s = Taylor(sin(a.taylor_vec[1]),gr)
    c = Taylor(cos(a.taylor_vec[1]),gr)
    # Promovemos los arreglos a y s,c.
    a, s = promote(a,s)
    a, c = promote(a,c)
    # El siguiente for es la relaci�n de recurrencia hallada que
    # usa los valores S_k,C_k previamente obtenidos.
    for k in range(1,gr)
        s.taylor_vec[k+1] = (1/k)*sum([c.taylor_vec[i]*(k-i+1)*a.taylor_vec[k+2-i] for i in range(1,k)])
        c.taylor_vec[k+1] = (-1/k)*sum([s.taylor_vec[i]*(k-i+1)*a.taylor_vec[k+2-i] for i in range(1,k)])
    end
    # Se regresa la estructura de Taylor con coeficientes calculados.
    return s
end

function cos(a::Taylor,gr::Int)
    # Este primer if est� armado para que no se regrese un error de
    # inexactitud de meter un entero a una exponencial que trabaja
    # m�nimo con n�meros flotantes.
    if eltype(a.taylor_vec) <: Integer
        a = Taylor(convert(Array{Float64,1}, a.taylor_vec))
    end
    # Inicializamos los Taylors s,cde grado gr.
    s = Taylor(sin(a.taylor_vec[1]),gr)
    c = Taylor(cos(a.taylor_vec[1]),gr)
    # Promovemos los arreglos a y s,c.
    a, s = promote(a,s)
    a, c = promote(a,c)
    # El siguiente for es la relaci�n de recurrencia hallada que
    # usa los valores S_k,C_k previamente obtenidos.
    for k in range(1,gr)
        s.taylor_vec[k+1] = (1/k)*sum([c.taylor_vec[i]*(k-i+1)*a.taylor_vec[k+2-i] for i in range(1,k)])
        c.taylor_vec[k+1] = (-1/k)*sum([s.taylor_vec[i]*(k-i+1)*a.taylor_vec[k+2-i] for i in range(1,k)])
    end
    # Se regresa la estructura de Taylor con coeficientes calculados.
    return c
end

# Implementamos la funci�n sin poner el orden expl�cito.
sin(a::Taylor) = sin(a,a.gr)
cos(a::Taylor) = cos(a,a.gr)


end
