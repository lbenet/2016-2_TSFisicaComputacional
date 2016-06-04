
#=
Éste archivo contiene los tests para el módulo *Taylor.jl*
Autores:
David Leonardo Galicia Praskauer
Arturo Cid Prati
Luis Yves Villegas Aguilar

12 de Abril de 2016
=#

include("Taylor.jl")
using Base.Test
using ADT

@test_throws MethodError (MethodError,Taylor("hi")) #Taylor sólo tienen Arrays como argumento, que pueden ser de n dims
@test typeof(Taylor([1,2]).pol)==Array{Int64,1}

@test Taylor([1.0,2.0])==Taylor([1,2,0]) #igualdad
@test Taylor([2])==Taylor(2) #taylor de un número
@test gradomax(Taylor([1,1,1]))==3 #comprobación del grado máximo del pol de taylor (sólo es la longitud del array)
@test gradomax(Taylor(8))==1 #un polinomio de un grado constante da un array de long. 1 (pol. de grado 0)

#promoción
@test prom(Taylor([1,2]),Taylor([1,2,3]))==Taylor([1,2,0])
@test prom(Taylor([1,2]),5)==Taylor([1,2,0,0,0])

#evaluación
@test evaluar(Taylor([1,1]),2)==3
@test evaluar(Taylor([1,2]),7)==evaluar(prom(Taylor([1,2]),20),7)

#suma
@test Taylor([1,2])+Taylor([1,2,0])==Taylor([2,4,0])
@test Taylor([1,2])+Taylor([3,4])==Taylor([3,4])+Taylor([1,2])
@test Taylor([1,2])+π==Taylor([1+π,2])

#resta
@test Taylor([1,2])-Taylor([1,2,0])==Taylor([0,0,0])
@test Taylor([1,2])-Taylor([3,4])==-(Taylor([3,4])-Taylor([1,2]))
@test Taylor([1,2])-π==Taylor([1-π,2])
@test Taylor([1,2])-π==-(π-Taylor([1,2]))

#multiplicación
@test Taylor([1,1])*Taylor([1,1])==Taylor([1,2,1])
@test π*Taylor([2,3])==Taylor([2*π,3*π])
@test π*Taylor([2,3])==Taylor([2,3])*π

#división
@test Taylor([6,3])/3==Taylor([2,1])
@test Taylor([1])==1/Taylor([1])
@test 1/prom(Taylor([1,-1]),10)==Taylor([1,1,1,1,1,1,1,1,1,1]) #serie geométrica hasta el 10o término
@test Taylor([0,0,1])/Taylor([0,1])==Taylor([0,1,0]) #x^3/x^2=x
@test Taylor([0,0,0,0,1])/Taylor([0,0,0,18])==Taylor([0,1])*(1/18)
@test 1/prom(Taylor([1,0,1]),5)==Taylor([1.0,0.0,-1.0,0.0,1.0]) #serie de 1/(1+x^2)

#exponencial
@test exp(Taylor([1]))==Taylor([e*1.0])
@test exp(Taylor(-Inf))==Taylor([0])
@test exp(Taylor([-1]))==Taylor(1/e) 
@test exp(prom(Taylor([0,1]),5))==Taylor([1,1,1/factorial(2),1/factorial(3),1/factorial(4)])
     #serie de taylor de exp(x) de orden 4 (5 términos)
#@test exp(prom(Taylor([0,2]),5))==exp(prom(Taylor([0,1]),5))*exp(prom(Taylor([0,1]),5))

#logaritmo natural
@test log(Taylor([1]))==Taylor([0])
@test log(Taylor([e]))==Taylor([1])
@test log(Taylor([0]))==Taylor([-Inf])

@test exp(log(Taylor([1,2])))==Taylor([1,2,0]) 
@test log(exp(Taylor([1,2])))==Taylor([1,2,0])
@test log(prom(Taylor([1,1]),7))==Taylor([0,1,-1/2,1/3,-1/4,1/5,-1/6]) #serie de taylor de log(1+x) con 7 términos

#potencias
@test Taylor([2,3,4])^1==Taylor([2,3,4])
@test Taylor([1,1])^6==Taylor([1,6,15,20,15,6,1]) #binomio de Newton (1+x)^6
@test prom(Taylor([1,1]),5)^(1/2)==Taylor([1,1//2,-1//8,1//16,-5//128]) #serie de Taylor de sqrt(1+x)
@test prom(Taylor([1,1]),4)^(-1/3) ==Taylor([1,-1/3,2/9,-14/81 #serie de 1/cubrt(1+x)


