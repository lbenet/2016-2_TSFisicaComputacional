__precompile__(true)

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
