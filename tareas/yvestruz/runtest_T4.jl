
#=
Este archivo contiene los tests para el módulo *AutomDiff*
Autores:
Luis Yves Villegas Aguilar
David Leonardo Galicia Praskauer
Arturo Cid Prati

15 de Marzo de 2016
=#

include("AutomDiff.jl")
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


