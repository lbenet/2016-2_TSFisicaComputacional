#Ejercicio 9.6.4

using ChaosCrypto
using PyPlot

function lorenz21(xx,t)
    x, y, z, xr, yr, zr= xx
    # x,y,z gobernados por: ρ = 60, σ = 10, β = 8/3 (default)
    # xr,yr,zr gobernados por parametros de ejercicio 9.6.2
    # Comenzamos a integrar numericamente
    
    AD.Taylor(x, [LO.σ])*(AD.Taylor(x,[y])-AD.Taylor(x,[x])),
    AD.Taylor(y,[x])*(AD.Taylor(y,[LO.ρ])-AD.Taylor(y,[z]))-AD.Taylor(y,[y]),
    AD.Taylor(z,[x])*AD.Taylor(z,[y])-AD.Taylor(z,[LO.β])*AD.Taylor(z,[z]),
    AD.Taylor(xr,[LO.σ])*(AD.Taylor(xr,[yr])-AD.Taylor(xr,[xr])),
    AD.Taylor(yr,[LO.ρ])*AD.Taylor(yr,[xr])-AD.Taylor(yr,[yr])-AD.Taylor(yr,[xr])*AD.Taylor(yr,[z]),
    AD.Taylor(zr,[x])*AD.Taylor(zr,[y])-AD.Taylor(zr,[LO.β])*AD.Taylor(zr,[z])
end

# Escogiendo condiciones iniciales diferentes para x, xr; y, yr..
vecs21,t21 = LO.integrador([1.,1.,1.,100,100,1.],lorenz21,10.)

function lorenz22(xx,t)
    x, y, z, xr, yr, zr = xx
    
    AD.Taylor(x, [LO.σ])*(AD.Taylor(x,[y])-AD.Taylor(x,[x])),
    AD.Taylor(y,[x])*(AD.Taylor(y,[LO.ρ])-AD.Taylor(y,[z]))-AD.Taylor(y,[y]),
    AD.Taylor(z,[x])*AD.Taylor(z,[y])-AD.Taylor(z,[LO.β])*AD.Taylor(z,[z]),
    AD.Taylor(xr,[LO.σ])*(AD.Taylor(xr,[yr])-AD.Taylor(xr,[xr])),
    AD.Taylor(yr,[LO.ρ])*AD.Taylor(yr,[x])-AD.Taylor(yr,[y])-AD.Taylor(yr,[x])*AD.Taylor(yr,[z]),
    AD.Taylor(zr,[xr])*AD.Taylor(zr,[yr])-AD.Taylor(zr,[LO.β])*AD.Taylor(zr,[zr])
    
end

vecs22,t22 = LO.integrador([1., 1., 1., 100, 1., 100],lorenz22,10.)

include("testejercicio2.jl")
title("Ejercicio 9.6.4 - 1, X vs Xr")
plot(t21, [x[1] for x in vecs21], label="x")
plot(t21, [x[4] for x in vecs21], label="xr")
fig = figure(1)

include("testejercicio2.jl")
title("Ejercicio 9.6.4 - 1, Y vs Yr")
plot(t21, [x[2] for x in vecs21], label="y")
plot(t21, [x[5] for x in vecs21], label="yr")
fig = figure(2)