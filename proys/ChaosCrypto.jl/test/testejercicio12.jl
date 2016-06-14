include("testejercicio1.jl")

# Graficamos la proyeccion (y,z) de ambas trayectorias
title("Ejercicio 9.6.3 - 2")
plot([x[2] for x in vecs],[x[3] for x in vecs])
xlabel("y")
ylabel("z")
legend()
show()

# title("Ejercicio 1.2")
# plot3D([x[1] for x in vecs],[x[2] for x in vecs],[x[3] for x in vecs])
# legend()
# show()