# El modulo en cuestion implementa diferenciacion automatica en Julia,
# definiendo a los duales y a sus operaciones
# Autores: Yuriko Yamamoto, Ignacio Vargas
# Fecha: 3 de mayo, 2016

# La siguiente instrucción sirve para *precompilar* el módulo
__precompile__(true)

module AD
    import Base: exp, log, sin, cos, tan, cot, sec, csc, asin, acos, atan, acot, asec, acsc, sinh, cosh, tanh, coth, sech, csch, asinh, acosh, atanh, acoth, asech, acsc


    exp(a::Dual) = Dual(exp(a.fun), a.der * exp(a.fun)) 

    log(a::Dual) = Dual(log(10,a.fun), a.der / a.fun) 

    sin(a::Dual) = Dual(sin(a.fun), a.der * cos(a.fun)) 

    cos(a::Dual) = Dual(cos(a.fun), -a.der * sin(a.fun))

    tan(a::Dual) = Dual(tan(a.fun), a.der * sec^2(a.fun))

    cot(a::Dual) = Dual(cot(a.fun), -a.der * csc^2(a.fun))

    sec(a::Dual) = Dual(sec(a.fun), a.der * sec(a.der) * tan(a.der))

    csc(a::Dual) = Dual(csc(a.fun), -a.der * csc(a.der) * cot(a.der))

    asin(a::Dual) = Dual(asin(a.fun), a.der * 1/sqrt(1 - (a.der)^2))

    acos(a::Dual) = Dual(acos(a.fun), -a.der * 1/sqrt(1 - (a.der)^2))

    atan(a::Dual) = Dual(atan(a.fun), a.der * 1/(1 + (a.der)^2))

    acot(a::Dual) = Dual(acot(a.fun), -a.der * 1/(1 + (a.der)^2))

    asec(a::Dual) = Dual(asec(a.fun), a.der * 1/(abs(a.der) * sqrt((a.der)^2 -1)))

    acsc(a::Dual) = Dual(acsc(a.fun), -a.der * 1/(abs(a.der) * sqrt((a.der)^2 -1)))
 
    sinh(a::Dual) = Dual(sinh(a.fun), a.der * cosh(a.fun))
 
    cosh(a::Dual) = Dual(cosh(a.fun), a.der * sinh(a.fun))

    tanh(a::Dual) = Dual(tanh(a.fun), a.der - a.der*tanh^2(a.fun))

    coth(a::Dual) = Dual(coth(a.fun), -a.der * (csch(a.der))^2)

    sech(a::Dual) = Dual(sech(a.fun), -a.der * sech(a.der) * tanh(a.der))

    csch(a::Dual) = Dual(csch(a.fun), -a.der * csch(a.der) *coth(a.der))

    asinh(a::Dual) = Dual(asinh(a.fun), a.der * 1/sqrt(1 + (a.der)^2))

    acosh(a::Dual) = Dual(acosh(a.fun), a.der * 1/sqrt((a.der)^2 -1))
    acosh(a::Dual) = Dual(acosh(a.fun), -a.der * 1/sqrt((a.der)^2 -1))

    atanh(a::Dual) = Dual(atanh(a.fun), a.der * 1/(1 - (a.der)^2))
    atanh(a::Dual) = Dual(atanh(a.fun), -a.der * 1/(1 - (a.der)^2))

    acoth(a::Dual) = Dual(acoth(a.fun), a.der * 1/(1-(a.der)^2))

    asech(a::Dual) = Dual(asech(a.fun), -a.der * 1/(abs(a.der) * sqrt(1 - (a.der)^2))) 
    asech(a::Dual) = Dual(asech(a.fun), a.der * 1/(abs(a.der) * sqrt(1 - (a.der)^2))) 

    acsch(a::Dual) = Dual(acsch(a.fun), -a.der * 1/(abs(a.der) * sqrt(1 + (a.der)^2)))

end
   