
#=
Este archivo contiene los tests para el módulo *AutomDiff*
Autores:
Luis Yves Villegas Aguilar
David Leonardo Galicia Praskauer
Arturo Cid Prati

15 de Marzo de 2016
=#

include("AutomDiff_T5.jl")
using Base.Test
using AD

#Consistencia en los tipos con la función xdual
@test typeof(xdual(1))==typeof(Dual(1,2)) #xdual crea un dual del mismo tipo 
@test xdual(2//1).fun ==xdual(2).fun #la comparación de las componentes de un dual se realiza correctamente
@test xdual(2//1).der ==xdual(2).der
@test typeof(xdual(1.2).der)==Float64
@test typeof(xdual(2).der)==Int64

#Tests para las operaciones con duales

@test Dual(1,2)==Dual(1.0,2.0) #duales se comparan entrada-entrada
@test Dual(1,2)+Dual(2,3)==Dual(3,5) #la suma es entrada-entrada
@test Dual(1,2)+Dual(2,3)==Dual(2,3)+Dual(1,2) #la suma es conmutativa también en duales
@test Dual(1,2)-Dual(2,3)==-(Dual(2,3)-Dual(1,2)) #duales respetan reglas de signos de reales
@test Dual(1,2)-Dual(2,3)==Dual(-1,-1) #la resta es entrada-entrada
@test Dual(1,2)==-Dual(-1,-2) #se pueden "factorizar" signos con duales

@test π*Dual(1,2)==Dual(π,2π)#distributividad
@test π*Dual(1,2)==Dual(1,2)*π #conmutatividad
@test π/Dual(1,2)==Dual(π/1,π/2) #real/dual
@test Dual(1,2)/π==Dual(1/π,2/π) #dual/real

f(x)=x^2
g(x)=x
@test f(Dual(3,1))==Dual(f(3),2*g(3)) #función aplicada a un dual
@test f(Dual(3,1))==g(Dual(3,1))*g(Dual(3,1))#derivada del producto
@test f(Dual(3,1))/g(Dual(3,1))==Dual(f(3)/g(3),(2*g(3)*g(3)-f(3))/(g(3))^2) #derivada del cociente
@test 1/f(Dual(0,1))==Dual(Inf,Inf) #1/0
@test 1/g(Dual(0,1))==Dual(Inf,1) 
#@test_throws Dual(3,2)===Dual(3,2)

#Tests para las definiciones de las funciones

#raíz cuadrada y cúbica
@test sqrt(xdual(2))==sqrt(Dual(2,1))
@test sqrt(xdual(2)).fun==sqrt(2)
@test sqrt(xdual(2)).der==1/(2*sqrt(2))
@test_throws(DomainError,sqrt(xdual(-1)))

@test cbrt(xdual(2))==cbrt(Dual(2,1))
@test cbrt(xdual(2)).fun==cbrt(2)
@test cbrt(xdual(2)).der==1/(3*cbrt(2)^2)
@test cbrt(xdual(-2))==Dual(cbrt(-2),1/(3*cbrt(-2)^2))
@test cbrt(xdual(-Inf))==Dual(-Inf,0.0)

@test sqrt(xdual(2))*cbrt(xdual(2))==Dual(sqrt(2)*cbrt(2),5/(6*(2^(1//6))))#derivada de x^1/2*x^1/3= derivada de x^1/6

#exponencial y logaritmo
@test exp(xdual(1))==Dual(e*1.0,e*1.0)
@test exp(xdual(π)).der==exp(xdual(π)).fun
@test exp(xdual(-Inf))==Dual(0.0,0.0)
@test 1/exp(xdual(1))==exp(Dual(-1,1))

@test log(xdual(1))==Dual(0,1)
@test log(xdual(0))==Dual(-Inf,Inf)
@test 5*log(xdual(5)).der==1
@test_throws(DomainError,log(xdual(-1)))

@test exp(log(xdual(π)))==Dual(π,1)
@test log(exp(xdual(π)))==Dual(π,1)


#funciones trigonométricas
@test cos(xdual(18))==Dual(cos(18),-sin(18))
@test cos(xdual(20))^3==Dual(cos(20)^3,-3cos(20)^2*sin(20))

@test sin(xdual(18))==Dual(sin(18),cos(18))
@test sin(xdual(20))^3==Dual(sin(20)^3,3sin(20)^2*cos(20))

@test tan(xdual(20))==Dual(tan(20),sec(20)^2)
@test cot(xdual(20))==Dual(cot(20),-csc(20)^2)

@test sec(xdual(20))==Dual(sec(20),tan(20)*sec(20))
@test csc(xdual(20))==Dual(csc(20),-cot(20)*csc(20))

#identidades pitagóricas
@test cos(xdual(20))^2+sin(xdual(20))^2==Dual(1,0)
@test tan(xdual(BigFloat(20.0)))^2+Dual(BigFloat(1),BigFloat(0))-sec(xdual(BigFloat(20.0)))^2==Dual(0,0)
@test cot(xdual(BigFloat(1)))^2+1-csc(xdual(BigFloat(1)))^2==Dual(0,0)

#funciones hiperbólicas
@test cosh(xdual(18))==Dual(cosh(18),sinh(18))
@test cosh(xdual(20))^3==Dual(cosh(20)^3,3cosh(20)^2*sinh(20))

@test sinh(xdual(18))==Dual(sinh(18),cosh(18))
@test sinh(xdual(20))^3==Dual(sinh(20)^3,3sinh(20)^2*cosh(20))

@test tanh(xdual(20))==Dual(tanh(20),sech(20)^2)
@test coth(xdual(20))==Dual(coth(20),-csch(20)^2)

@test sech(xdual(20))==Dual(sech(20),-tanh(20)*sech(20))
@test csch(xdual(20))==Dual(csch(20),-coth(20)*csch(20))

