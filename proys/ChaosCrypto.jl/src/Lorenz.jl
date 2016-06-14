__precompile__(true)
"""
#Modulo LO (Lorenz.jl)

``dx/dt = σ(y-x),``
``dy/dt = ρx - y - xz,``
``dz/dt = xy -bz.``

Operaciones que conciernen a las ecuaciones de Lorenz.
Los parámetros ρ, σ, b deben de ser reales y positivos.

El módulo trabaja por medio de integración con polinomios de Taylor.
"""

module LO
    export generarTaylor,generarSerie
    export generaIntervalo,horner,integrador
    export σ,ρ,β


"""
#Generador Taylor

        generarTaylor(condIni, funcion, t)

Genera polinomio Taylor dadas condiciones iniciales, una función, y parámetro tiempo.
"""

    function generarTaylor(condIni, funcion, t)
        funcion(condIni, t)
    end

"""
#Generador Serie

    generarSerie(polTalor)

Genera un arreglo x de números flotantes. Pone los valores en los que se calcula la serie en x, y copia a arreglo valores
a los coeficientes. Regresa a x con los valores de los coeficientes entre su orden.

"""

    function generarSerie(polTalor)
        x = Float64[]
        push!(x, polTalor.ini)
        valores = copy(polTalor.coef)
        for i in 1:length(valores)
            push!(x, valores[i]/i)
        end
        x
    end


    function generaIntervalo(lista)
        p = length(lista)
        h = lista[end]
        while (lista[p] == 0)
            p -= 1
            h = lista[p]
        end
        ϵ = 1e-3
        (ϵ/h)^(1/p)
    end

"""
#Método Horner

    horner(x,h = 1e-3)

Simple implementación del método de Horner, con h = 1e-3 por default.

"""

    function horner(x,h = 1e-3)
        p = length(x)
        xt = zeros(x)
        xt[p-1] = x[p-1] + h*x[p]
        for i in 3:p
            xt[p-i+1] = x[p-i+1] + h*xt[p-i+1+1]
        end
        xt[1]
    end

    σ = 10
    ρ = 60
   #ρ = 28
    β = 8/3

"""
## Integrador por Taylor

        integrador(x0, f, tf)

Integración por método de Taylor. La utilizamos para resolver las ecuaciones de Lorenz.
Se vale de las funciones generarTaylor, generarSerie y horner - con parámetros de condiciones iniciales (x0), 
función, y tiempo final.

"""
    function integrador(x0, f, tf)
        a = generarTaylor(x0,f,0.)
        b = map(generarSerie,a)
        suma = map(horner,b)
        sol = Array{Float64,1}[x0, [suma...]]
        t = [0.,1e-3]
        while t[end] < tf
            a = generarTaylor(sol[end],f,t[end])
            b = map(generarSerie,a)
            suma = map(horner,b)
            push!(sol,[suma...])
            push!(t,t[end]+1e-3)
        end
        sol,t
    end
end
