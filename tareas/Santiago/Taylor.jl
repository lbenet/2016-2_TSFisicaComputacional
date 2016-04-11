__precompile__(true)

module TS
	import Base: +, -, *, /, ^, ==
	import Base: exp, log, sin, cos

	export Taylor
"""
Taylor

Es un tipo parametrizado por subtipos de Number. Consiste en una entrada
vecotrial Taylor.taylor_vec donde se mapean los primeros k coeficientes
de Taylor. La k escogida en esta estructura en particular es 10 y est
implementada para que si se mete un arreglo con menos de 10 entradas se
complemente con ceros restantes y si se mete un arreglo con ms de 10 
entradas se tomen slo las primeras 10.
"""
type Taylor{T<:Number}
    # Aqu se declara taylor_vec que es un Array unidimensional de tipo T.
    taylor_vec::Array{T, 1}
    
    # Aqu hay un constructor interno que lo que va a hacer es tomar el
    # arreglo taylor_vec y ver si tiene un tamao mayor o menor que la
    # k que es donde se trunca la serie (en este caso k = 10).
    function Taylor(taylor_vec::Array{T, 1})
        if size(taylor_vec, 1) < 10
            # Si hay menos de k entradas, se utiliza vcat para meter ceros 
            # restantes en taylor_vec y luego se reconstruye.
            taylor_vec = vcat(taylor_vec, zeros(T, 10 - size(taylor_vec,1)))
            new(taylor_vec)
            else
            # En caso contrario se seleccionan las primeras k entradas.
            taylor_vec = taylor_vec[1:10]
            new(taylor_vec)
        end
    end
end

# Este es un constructor externo para evitar problemas de conversin. Vase:
# http://docs.julialang.org/en/release-0.4/manual/constructors/
# Especficamente la seccin de Parametric Constructors.
Taylor{T<:Number}(taylor_vec::Array{T, 1}) = Taylor{T}(taylor_vec)

# Adicionalmente se declara un mtodo para tener Taylor para un escalar.
Taylor{T<:Number}(a::T) = Taylor(a*ones(T,1))

# Ahora se definen las operaciones artmeticas:
# La suma entre los dos Taylors y entre un escalar.
+(a::Taylor, b::Taylor) = Taylor(a.taylor_vec + b.taylor_vec)
+(a::Taylor, b::Number) = a + Taylor(b)
+(a::Number, b::Taylor) = Taylor(a) + b
# La resta se implementa igualmente.
-(a::Taylor, b::Taylor) = Taylor(a.taylor_vec - b.taylor_vec)
-(a::Taylor, b::Number) = a - Taylor(b)
-(a::Number, b::Taylor) = Taylor(a) - b
# Los operadores unitarios y la comparacin.
+(a::Taylor) = a
-(a::Taylor) = Taylor(-a.taylor_vec)
==(a::Taylor, b::Taylor) = a.taylor_vec == b.taylor_vec
# La mutiplicacin se implementa usando dos for's. El primero que va
# desde 1, size(a.taylor_vec, 1) construye el arreglo entero de k 
# elementos. El segundo for anidado es la suma que se da en la teora
# de arriba. Slo hay que tener cuidado porque los ndices empiezan en
# uno y no en cero como en la teora.
*(a::Taylor, b::Taylor) = Taylor([sum([a.taylor_vec[i]*b.taylor_vec[k-i+1] 
        for i in range(1,k)]) for k in range(1,size(a.taylor_vec,1))])
*(a::Taylor, b::Number) = a*Taylor(b)
*(a::Number, b::Taylor) = Taylor(a)*b

# La divisin resulta ser un poco ms complicada as que se opta por
# definirla en una funcin explcita y luego meterla en el mtodo de /.
function /(a::Taylor, b::Taylor)
    # Lo primero que se hace es construir el arreglo de la divisin que
    # sea del mismo tamao de k elementos y lleno de ceros enteros.
    # Enteros para que se puedan promover en el caso de que los arreglos
    # sean de tipo Float, Complex, etc.
    div_vec = 0*a.taylor_vec
    # El primer elemento es sencillamente:
    div_vec[1] = a.taylor_vec[1]/b.taylor_vec[1]
    # De esta manera nos evitamos problemas con la suma que ira sobre un
    # rango 1 a 0.
    # Ahora podemos llenar los elementos que falta. Nuestro rango empieza
    # en 2 pero es de longitud(div_vec,1) -1:
    for k in range(2, size(div_vec,1)-1)
        # sd ser la suma que aparece en la ecuacin, inicializada como cero.
        sd = 0
        # Con el siguiente for se ejecuta la suma utilizando los coeficientes
        # de div_vec que ya se generaron.
        for i in range(1,k-1)
            sd = sd + (div_vec[i])*b.taylor_vec[k-i+1]
        end
        div_vec[k] = (1/b.taylor_vec[1])*(a.taylor_vec[k]-sd)
    end
    return Taylor(div_vec)
end       
/(a::Taylor, b::Number) = a/Taylor(b)
/(a::Number, b::Taylor) = Taylor(a)/b

# Esta es la funcin exponencial aplicada a estructuras Taylor.
function exp(a::Taylor)
    # Este primer if est armado para que no se regrese un error de
    # inexactitud de meter un entero a una exponencial que trabaja 
    # mnimo con nmeros flotantes.
    if eltype(a.taylor_vec) <: Integer
        a = Taylor(convert(Array{Float64,1}, a.taylor_vec))
    end
    # Inicializamos el arreglo e_vec que ser quien llene la 
    # estructura Taylor asociada a la exponencial:
    e_vec = 0.*a.taylor_vec
    # Ahora la condicin inicial sobre el arreglo.
    e_vec[1] = exp(a.taylor_vec[1])
    # El siguiente for es la relacin de recurrencia hallada que
    # usa los valores E_k previamente obtenidos.
    for k in range(2,size(e_vec,1)-1)
        e_vec[k] = (1/(k-1))*sum([(k-j)*e_vec[j]*a.taylor_vec[k-j+1] for j in range(1,k)])
    end
    # Se regresa la estructura de Taylor con coeficientes calculados.
    return Taylor(e_vec)
end

# Esta es la funcin logaritmo aplicada a estructuras Taylor.
function log(a::Taylor)
    # Este primer if est armado para que no se regrese un error de
    # inexactitud de meter un entero a un logaritmo que trabaja 
    # mnimo con nmeros flotantes.
    if eltype(a.taylor_vec) <: Integer
        a = Taylor(convert(Array{Float64,1}, a.taylor_vec))
    end
    # Inicializamos el arreglo e_vec que ser quien llene la 
    # estructura Taylor asociada al logaritmo:
    l_vec = 0.*a.taylor_vec
    # Ahora la condicin inicial sobre el arreglo.
    l_vec[1] = log(a.taylor_vec[1])
    # Ahora vamos a crear el arreglo de la derivada de a.
    b = Taylor(convert(Array{eltype(l_vec),1},[i*a.taylor_vec[i+1] for i in range(1,size(a.taylor_vec,1)-1)]))
    # El siguiente for es la relacin de recurrencia hallada que
    # es bsicamente el de la divisin de dos series de Taylor.
    for k in range(1,size(l_vec,1)-1)
        l_vec[k+1] = (1/k)*(b/a).taylor_vec[k]
    end
    # Se regresa la estructura de Taylor con coeficientes calculados.
    return Taylor(l_vec)
end

# Con el logaritmo y la exponencial es muy sencillo poner la potencia:
# Recuerdese que Julia pide primero declararla para exponentes enteros:
function ^(a::Taylor,b::Integer)
    c = a
    for i in range(1,b-1)
        c = c*a
    end
    return c
end
^(a::Taylor,b::Number) = exp(b*log(a))

# Teniendo la potencia a nuestra disposicin es casi trivial escribir
# las funciones coseno y seno en trminos de su serie de pontecias.
cos(a::Taylor) = sum([((-1)^(k)/(factorial(2*k)))*a^(2*k) for k in range(0,size(a.taylor_vec,1)-1)])
sin(a::Taylor) = sum([((-1)^(k)/(factorial(2*k+1)))*a^(2*k+1) for k in range(0,size(a.taylor_vec,1)-1)])
end
