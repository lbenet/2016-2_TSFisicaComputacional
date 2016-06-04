using WAV
using PyPlot
using ChaosCrypto

# cosas para el audio
sonido,fs = wavread("test.wav");
valores = [x[1] for x in sonido]
usados = copy(valores);

function lorenzAudio(xx,t)
    x, y, z, xr, yr, zr = xx
    
    m = shift!(valores)
    s = x + m
        [AD.Taylor(x,[LO.σ*(y-x)]), 
        AD.Taylor(y,[(LO.ρ*x - y - x*z)]), 
        AD.Taylor(z,[x*y-LO.β*z]),
        AD.Taylor(xr,[LO.σ*(yr-xr)]), 
        AD.Taylor(yr,[(LO.ρ*s - yr - s*zr)]), 
        AD.Taylor(zr,[(s*yr - LO.β*zr)])]
end

vecs,t = LO.integrador([1.,1.,1.,1.,100.,100.],lorenzAudio,(length(usados)-2)*1e-3)

m=usados[1:length(vecs)]
s=[x[1] for x in vecs]+m
mhat = s-[x[4] for x in vecs]
plot(t,m,".")
plot(t,mhat,".")


plot(mhat)
ylim(-.2,.3)


plot(usados)
