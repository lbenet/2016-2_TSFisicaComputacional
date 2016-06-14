include("testejercicio1.jl")

# Graficamos a y(t) y yr(t) en el mismo plot. Observamos que eventualmente las dos series de tiempo se mergen.
title("Ejercicio 9.6.3 - 1")
plot(t, [x[2] for x in vecs], label="y")
plot(t, [x[5] for x in vecs], label="yr")
#plot(t, [x[2] for x in vecs], ".", label=L"$y$")
#plot(t, [x[5] for x in vecs], label=L"$y_r$")
xlabel("t")
ylabel("y, yr")
legend()
show()