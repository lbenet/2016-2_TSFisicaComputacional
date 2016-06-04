#Ejercicio 9.6.3

using ChaosCrypto
using PyPlot

function lorenz1(xx,t)
    x, y, z, xr, yr, zr= xx
    # x,y,z gobernados por: ρ = 60, σ = 10, β = 8/3 (default)
    # xr,yr,zr gobernados por parametros de ejercicio 9.6.2
    # Comenzamos a integrar numericamente
    AD.Taylor(x, [LO.σ])*(AD.Taylor(x,[y])-AD.Taylor(x,[x])),
    AD.Taylor(y,[x])*(AD.Taylor(y,[LO.ρ])-AD.Taylor(y,[z]))-AD.Taylor(y,[y]),
    AD.Taylor(z,[x])*AD.Taylor(z,[y])-AD.Taylor(z,[LO.β])*AD.Taylor(z,[z]),
    AD.Taylor(xr,[LO.σ])*(AD.Taylor(xr,[y])-AD.Taylor(xr,[xr])),
    AD.Taylor(yr,[LO.ρ])*AD.Taylor(yr,[xr])-AD.Taylor(yr,[yr])-AD.Taylor(yr,[xr])*AD.Taylor(yr,[zr]),
    AD.Taylor(zr,[xr])*AD.Taylor(zr,[yr])-AD.Taylor(zr,[LO.β])*AD.Taylor(zr,[zr])
end

# Escogiendo condiciones iniciales diferentes para y, yr; z, zr..
vecs,t = LO.integrador([1.,1.,1.,1.,100.,100.],lorenz1,10.);