using ChaosCrypto
using PyPlot
using Interact
σ=16.
ρ=45.6
β=4.
fig = figure()
@manipulate for φ in .001:.01:.1
    withfig(fig) do

function lorenzSin(xx,t)
    x, y, z, xr, yr, zr = xx
    
    m = sin(φ*t)
    s = x + m
    
    [σ*(y-x), (ρ*x - y - x*z), x*y-β*z,
     σ*(yr-xr), (ρ*s - yr - s*zr) , (s*yr - β*zr)]
end

xs, ts = RK.integrar(lorenzSin,[1., 1., 1., 10., 10., 10.],0. ,100. ,1e-3);
m=sin(ts)
s=[x[1] for x in xs]+m
mhat = s-[x[4] for x in xs]
title("Ejercicio 3 con RK")
plot(ts,m)
plot(ts,mhat)
end
end
