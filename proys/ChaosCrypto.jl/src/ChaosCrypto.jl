__precompile__(true)

"""
#Módulo ChaosCrypto

Módulo principal que llama al código de los módulos de Diferenciación Automática, Runge-Kutta, Integración por Taylor
y el despliegue de audio en los Jupyter notebooks. Busca simplificar la llamada de estos módulos juntándolos todos en uno.

Véanse

```Automata.jl
Lorenz.jl
RK4.jl
AudioDisplay.jl
```

"""

module ChaosCrypto

export AD, LO, RK, audioplayer

export Taylor, paso2,paso1
export igualdad,logo,expo,seno,coseno

export generarTaylor,generarTaylor,generarSerie
export generaIntervalo,horner,integrador
export σ,ρ,β

export runge4, integrar

include("Automata.jl")
include("Lorenz.jl")
include("RK4.jl")
include("AudioDisplay.jl")

end
