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
function gradomax(a::Taylor,b::Taylor)
    return max(length(a.pol),length(b.pol))
end
gradomax(a::Taylor) = length(a.pol)

# Promoción
function prom(a::Taylor,b::Taylor = a)
    q = [a.pol; fill(0,gradomax(a,b)-gradomax(a))]
    return Taylor(q)
end
prom(a::Taylor,n::Integer) = prom(a,Taylor(zeros(n)))

#Evaluación
"""Función que evalúa un Taylor en el punto x0 (debe ser cercano al punto alrededor del que se construye el polinomio de Taylor)
"""
function evaluar(a::Taylor,x0)
    n = gradomax(a);
    ex = :(0)
    
    for k = 1:n
        ex = :($ex + $a.pol[$k]*$x0^$(k-1))
    end
    return eval(ex)
end

import Base: +, -, *, ^, /, ==
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
function /(b::Real, a::Taylor)
    n = gradomax(a);
    r = Taylor(zeros(n));
    A = prom(a,r);
    s = 1; # índice desde donde empezamos
    if A.pol[s] == 0 # checamos si el primer término no es nulo
        s += 1;
    end
    r.pol[1] = b/A.pol[s];
    for k = s:n-1
        suma = 0;
        for j = 0:k-1
            suma += r.pol[j+1]*A.pol[k-j+1]
        end
        r.pol[k+1] = (-suma)/A.pol[s];
    end
    return r
end
/(a::Taylor, b::Taylor) = a*(1/b)
/(a::Taylor, k::Number) = Taylor(a.pol/k)
/(k::Number, a::Taylor) = Taylor(k)/a

# Igualdad
==(a::Taylor, b::Taylor) = a.pol == b.pol

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
    while a.pol[s] == 0
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

# Potencia
function ^(a::Taylor, n::Integer)
    ex = :($a)
    k = 1;
    while k < n
        ex = :($ex * $a)
        k += 1
    end
    return eval(ex)
end
^(a::Taylor, n::Number) = Taylor(exp(n*log(a)))

# Seno

# Coseno
