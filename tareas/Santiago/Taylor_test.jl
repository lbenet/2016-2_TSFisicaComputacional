# Este archivo tiene los tests del modulo TS (Taylor Series)

using Base.Test
using TS

# Primero vamos a probar que nuestra estructura Taylor tiene solo espacio para 10 entradas. Esto se puede cambiar
# muy rapido en el codigo cambiando solo la definicion de Taylor.

@test Taylor([1,1,1]) == Taylor([1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,2,1])

# Ahora corresponden las pruebas de operaciones aritmetica basica.
@test Taylor([1,0]) + Taylor([0,1]) == Taylor([1,1])
@test Taylor([0,1]) + 1 == Taylor([1,1])
@test Taylor(ones(8))/Taylor([1,1]) * Taylor([1,1]) == Taylor(ones(8))
@test Taylor(ones(5)) == (-1)*(-Taylor(ones(5)))
@test Taylor([1.,2.,1.])/Taylor([1.,1.]) == Taylor([1.,1.])

# Despues podemos probar ya las funciones implentadas.
a = Taylor([1,2,1])
@test exp(a).taylor_vec[2] == 2*e
@test exp(a).taylor_vec[3] == 3*e
@test log(Taylor(1)) == Taylor([0])
@test log(a).taylor_vec[3] == -1.
@test exp(log(a)).taylor_vec[1] == 1.0
@test a^2 == a*a
@test (a^(1/2)).taylor_vec[1:2] == [1.,1.]
@test (a^(2im)).taylor_vec[3] == -(8.0+2.0im)
@test cos(a).taylor_vec[1] == cos(1)
@test ((cos(a))^2 + (sin(a))^2).taylor_vec[1] == 1.0
