#Ejercicio 9.6.5

using ChaosCrypto
using PyPlot

fig = figure()

@manipulate for ρr in [1,10,24,25,60]
    withfig(fig) do
        function lorenzSin(xx,t)
            x, y, z, xr, yr, zr = xx
            
            φ = 1e-3
            m = sin(φ*t)
            s = x + m
            
            [AD.Taylor(x,[LO.σ*(y-x)]), 
                AD.Taylor(y,[(ρr*x - y - x*z)]), 
                AD.Taylor(z,[x*y-LO.β*z]),
                AD.Taylor(xr,[LO.σ*(yr-xr)]), 
                AD.Taylor(yr,[(ρr*s - yr - s*zr)]), 
                AD.Taylor(zr,[(s*yr - LO.β*zr)])]
        end
        vecs,t = LO.integrador([1.,1.,1.,1.,100.,100.],lorenzSin,10.)

        m = sin(t)
        s = [x[1] for x in vecs] + m
        mhat = s-[x[4] for x in vecs]
        title("Comparación señal y resultado")
        plot(t,m, label="Señal")
        plot(t,mhat, label="Resultado")
        xlabel("tiempo")
        ylabel("amplitud")
        legend()
    end
end