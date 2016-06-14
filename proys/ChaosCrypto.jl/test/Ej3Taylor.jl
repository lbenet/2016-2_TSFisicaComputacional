#Ejercicio 9.6.5
#Variando frecuencias

push!(LOAD_PATH, "/Users/Yuriko/ChaosCrypto.jl/test/")

using ChaosCrypto
using PyPlot
using Interact

export lorenzSinTaylor

fig = figure()
@manipulate for φ = .0001:.01:.1 
    withfig(fig) do
    @show(φ)



function lorenzSinTaylor(xx, t) 
    x, y, z, xr, yr, zr = xx #6 ecuaciones de Lorenz
    
   
    m = sin(φ*t) #Mensaje, señal
    s = x + m
    
        #Entrada 1 de xx, correspondiente a la ecuación x
        [AD.Taylor(x, [LO.σ *(y - x)]),
        #Entrada 2 de xx, correspondiente a la ecuación y
        AD.Taylor(y, [LO.ρ * x - y - x * z]),
        #Entrada 3 de xx, correspondiente a la ecuación z
        AD.Taylor(z, [x * y - LO.β * z]),
        #Entrada 4 de xx, correspondiente a la ecuación xr
        AD.Taylor(xr, [LO.σ * (yr - xr)]),
        #Entrada 5 de xx, correspondiente a la ecuación yr
        AD.Taylor(yr, [LO.ρ * s - yr - s * zr]),
        #Entrada 6 de xx, correspondiente a la ecuación zr
        AD.Taylor(zr, [s * yr - LO.β * zr])]

end

xs, ts = LO.integrador([1.0,1.0,1.0,1.0,10.0,10.0], lorenzSinTaylor, 500.0)

    m = sin(ts) #señal
    s = [x[1] for x in xs] + m #s = x + m
mhat = s - [x[4] for x in xs] #mensaje recibido, que es s - xr, debe ser aproximadamente m 
    

   



end #del manipulate
end