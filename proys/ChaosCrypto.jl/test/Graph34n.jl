include("Ej3Taylor.jl")

println(ts)

#Gráfica 1
plot(ts,m) #Gráfica de la señal o mensaje
plot(ts,mhat) #Gráfica de la señal recibida
title("Gráfica de la señal enviada y recibida")
xlabel("t")
ylim(-5, 5)
xlim(0.0, 200.0)

