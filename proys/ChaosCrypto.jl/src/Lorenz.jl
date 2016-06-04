#Modulo Lorenz
__precompile__(true)

module LO
    export generarTaylor,generarSerie
    export generaIntervalo,horner,integrador
    export σ,ρ,β

    function generarTaylor(condIni, funcion, t)
        funcion(condIni, t)
    end

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
