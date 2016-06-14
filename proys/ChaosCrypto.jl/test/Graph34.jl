#Ejercicio 9.6.5
#Variando frecuencias

push!(LOAD_PATH, "/Users/Yuriko/ChaosCrypto.jl/test/")


using ChaosCrypto
using PyPlot
using Interact


fig = figure()
@manipulate for φ = .0001:.005:.3 
    withfig(fig) do
    #@show(φ)



function lorenzSinTaylor(xx, t) 
    x, y, z, xr, yr, zr = xx #6 ecuaciones de Lorenz
    
    #φ = 1e-3 
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
#φ = 1e-3
    m = sin(ts) #señal
    s = [x[1] for x in xs] + m #s = x + m
mhat = s - [x[4] for x in xs] #mensaje recibido, que es s - xr, debe ser aproximadamente m 
    

#Gráfica 1
plot(ts,m, label="mensaje enviado") #Gráfica de la señal o mensaje
plot(ts,mhat, label="mensaje recibido") #Gráfica de la señal recibida
ylim(-5, 5)
xlim(0.0, 200.0)
title("Gráfica de la señal enviada y recibida")
xlabel("Tiempo")
ylabel("amplitud")
legend()
#show()



end
end #del manipulate
