include("testejercicio2.jl")
title("Ejercicio 9.6.4 - 2, Z vs Zr")
plot(t22, [x[3] for x in vecs22], label="z")
plot(t22, [x[6] for x in vecs22], label="zr")
xlabel("t")
ylabel("z, zr")
legend()
show()