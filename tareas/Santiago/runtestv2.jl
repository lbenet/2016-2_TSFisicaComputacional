# Este archivo incluye los tests del módulo AD

using Base.Test
using AD

# A continuación vienen los tests que implementaron y que deben ser suficientemente exhaustivos

# Primero se ve que no se puede meter un tipo que no sea un subtipo de Real para el dual:
@test_throws StackOverflowError Dual("santiago", 1...)

# Luego vemos que el dual de una constante sí sea lo que se espera:
@test Dual(pi) == Dual(pi, 0...)

# Lo mismo para la funcion identidad xdual:
@test xdual(pi) == Dual(pi,1...)

# Podemos añadir pruebas del tipo
@test xdual(2 + pi) == Dual(2 + pi, 1...)

# Ahora hacemos pruebas de las operaciones aritméticas.
@test xdual(2) + xdual(pi) == Dual(2+pi, 1+1...)
@test xdual(pi) + 2 == Dual(pi+ 2, 1...)
@test xdual(2) - xdual(pi) == -(xdual(pi) - xdual(2))
@test xdual(2)*xdual(pi) == Dual(2*pi, 2*1+pi*1...)
@test xdual(2)/Dual(pi) == Dual(2/pi, (1-(2/pi)*0)/pi...)
@test (xdual(2)*5)*(1/5) == xdual(2)
@test xdual(pi)^1 == xdual(pi)


# TEST PARA LA SEGUNDA PARTE.
@test cos(xdual(pi)).der == -sin(pi)
@test exp(xdual(1.0)^2).der == 2exp(1.0)
@test sinh(xdual(pi)).der == cosh(pi)
@test exp(log(xdual(1.0))).der == xdual(1.0).der
@test atan(xdual(1)).der == 1/2
@test (cos(xdual(2))^2).der+((sin(xdual(2)))^2).der ==0

