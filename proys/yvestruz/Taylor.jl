#=
Este módulo contiene las bases para la diferenciación automática en Julia. Permite calcular de manera exacta la derivada de órdenes superiores utilizando el álgebra de polinomios truncados de Taylor, introduciendo un nuevo tipo de estructura: Taylor.

Autores:
Luis Yves Villegas Aguilar
David Leonardo Galicia Praskauer
Arturo Cid Prati

12 de Abril de 2016
=#

__precompile__(true)

module ADT
    import Base: +, -, *, /, ^, ==

    export Taylor, gradomax, prom, evaluar




type Taylor{T<:Number}
    # código:
    pol :: Array{T}
end

# Consistencia
Taylor(a::Number) = Taylor([a])
Taylor(a::Complex) = Taylor ([a])

# Grado máximo
""""
Regresa la longitud del array más grande entre dos polinomios. grado del polinomio x = gradomax(x)-1
"""
function gradomax(a::Taylor,b::Taylor)
    return max(length(a.pol),length(b.pol))
end
gradomax(a::Taylor) = length(a.pol)

# Promoción
"""prom(a::Taylor,b::Taylor = a )
función que promueve el primer polinomio a al grado del segundo b si este es mayor. De no se así regresa a.
"""
function prom(a::Taylor,b::Taylor = a)
    q = [a.pol; fill(0,gradomax(a,b)-gradomax(a))]
    return Taylor(q)
end
prom(a::Taylor,n::Integer) = prom(a,Taylor(zeros(n)))

#evaluación
"""Función que evalúa un Taylor en el punto x0 (debe ser cercano al punto alrededor del que se construye el polinomio de Taylor)
"""
function evaluar(a::Taylor,x0::Number)
    n = gradomax(a);
    ex = :(0)

    for k = 1:n
        ex = :($ex + $a.pol[$k]*$x0^$(k-1))
    end
    return eval(ex)
end


import Base: +, -, *, ^, /, ==

# Igualdad
function ==(a::Taylor, b::Taylor)
    A=prom(a,b)
    B=prom(b,a)
    return A.pol == B.pol
end

# Suma
+(a::Taylor, b::Taylor) = Taylor(prom(a,b).pol+prom(b,a).pol)
+(a::Taylor, k::Number) = a + Taylor(k)
+(k::Number, a::Taylor) = a + Taylor(k)

# Resta
-(a::Taylor, b::Taylor) = Taylor(prom(a,b).pol-prom(b,a).pol)
-(a::Taylor, k::Number) = a - Taylor(k)
-(k::Number, a::Taylor) = Taylor(k) - a
-(a::Taylor) = Taylor(-a.pol)

# Multiplicación.
function *(a::Taylor, b::Taylor)
    n = gradomax(a)+ gradomax(b)-1;
    r = Taylor(zeros(n));
    A = prom(a,r);
    B = prom(b,r)
    for k = 0:n-1
        suma = 0;
        for j = 0:k
            suma += A.pol[j+1]*B.pol[k-j+1];
        end
        r.pol[k+1] = suma;
    end
    return r
end
*(a::Taylor, k::Number) = Taylor(k*a.pol)
*(k::Number, a::Taylor) = Taylor(k*a.pol)

# División
function /(A::Taylor, B::Taylor)
    a = prom(A,B);
    b = prom(B,A);

    n = gradomax(a);

    r = Taylor(zeros(n));
    s = 1; # índice desde donde empezamos

    while b.pol[s] == 0 # checamos si el primer término no es nulo
        s += 1;
    end

    r.pol[1] = a.pol[s]/b.pol[s];

    for k = (s+1):n
        suma = 0;

        for j = 0:k-1
            suma += r.pol[j+1]*b.pol[k-j]
        end

        r.pol[k-s+1] = (a.pol[k]-suma)/b.pol[s];
    end
    return r
end
/(a::Taylor, k::Number) = Taylor(a.pol/k);
/(k::Number, a::Taylor) = Taylor(k)/a
div_ex(a::Taylor, b::Taylor, n::Integer) = prom(a,n)/prom(b,n)
div_ex(a::Taylor, k::Number, n::Integer) = prom(a,n)/k
div_ex(k::Number, a::Taylor, n::Integer) = Taylor(k)/prom(a,n)


# # # Funciones

import Base: exp, log, sin, cos

# Exponencial
function exp(a::Taylor)
    n = gradomax(a); # grado máximo
    exp_t = Taylor(zeros(n)); # prealocación de memoria
    exp_t.pol[1] = exp(a.pol[1]); # definimos el primer elemento de la serie
    for k = 2:n
        suma = 0;
        for j = 1:k
            suma += (k-j)*a.pol[k-j+1]*exp_t.pol[j];
        end
        exp_t.pol[k] = suma*(1/(k-1));
    end
    return exp_t
end
exp(a::Taylor,n::Integer) = exp(prom(a,n))

# Logaritmo
function log(a::Taylor)
    n = gradomax(a); # grado máximo
    L = Taylor(zeros(n));
    L.pol[1] = log(a.pol[1]);
    s = 1; # índice desde donde empezamos

    #hallar el índice del primer término no nulo o hacer s=n
    while s<=n && a.pol[s] == 0
        s += 1
    end

    for k = (s+1):n
        suma = 0;
        for j = (s+1): k
            suma += (j-1)*L.pol[j]*a.pol[k-j+1];
        end
        L.pol[k] = (1/a.pol[s])*(a.pol[k]-suma/(k-s))
    end
    return L
end
log(a::Taylor, n::Integer) = log(prom(a,n))

# Potencia de un Número Entero
function ^(a::Taylor, n::Integer)
    if n != 0
        ex = :($a)
        k = 1;
        while k < n
            ex = :($ex * $a)
            k += 1
        end
        return eval(ex)
    else
        return Taylor(ones(1))
    end
end

#Potencia de un Número Cualquiera
^(a::Taylor, n::Number) = exp(n*log(a))

# Seno
import Base: sin,cos

function sin(a::Taylor)
    n = gradomax(a);
    S = Taylor(zeros(n));
    for k = 0:9
        S += (-1)^k * a^(2*k+1) /factorial(2*k + 1);
    end
    return S
end

# Coseno
function cos(a::Taylor)
    n = gradomax(a);
        C = Taylor(zeros(n))
        for k = 0:9
            C += (-1)^k * a^(2*k) / factorial(2*k);
        end
    return C
end


end
