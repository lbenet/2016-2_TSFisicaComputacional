
# Aqui viene una explicación de lo que se hace en el módulo, los autores y la fecha
#Este módulo implementa la definición del tipo dual, así como sus operaciones.



__precompile__(true)

module AD
    import Base: +, -, *, /, ^

    export Dual, xdual

type Dual{D <: Real}
    fun::D
    der::D
    
end
Dual(a,b)=Dual(promote(a,b)...)

Dual(a)=Dual(a,0)

function xdual(x)
    if typeof(x) <: Real
        fun=x
        der=1
        return Dual(fun,der)
    end   
end

import Base: +, -, *, /, ^

# Aqui se implementan los métodos necesarios para cada función
for fn1 = (:+, :-)
    println(fn1)
    ex = quote
        function ($fn1)(a::Dual, b::Dual)
            xx = ($fn1)(a.fun, b.fun)
            yy = ($fn1)(a.der, b.der)
            return Dual(xx, yy)
        end
    end
    println(ex)
    @eval $ex
end

for fn2 in (:*,), fn3 in (:+,)
    println(fn2)
    ex = quote
        function ($fn2)(a::Dual, b::Dual)
            xx = ($fn2)(a.fun, b.fun)
            yy = ($fn3)(($fn2)(a.fun, b.der),($fn2)(a.der, b.fun))
            return Dual(xx, yy)
        end
    end
    println(ex)
    @eval $ex
end


for fn4 = (:/,), fn5 = (:*,), fn6 = (:-,)
    println(fn4)
    ex = quote
        function ($fn4)(a::Dual, b::Dual)
            xx = ($fn4)(a.fun, b.fun)
            yy = ($fn4)(($fn6)(a.der, ($fn5)(($fn4)(a.fun,b.fun),b.der)),b.fun)
            return Dual(xx, yy)
        end
    end
    println(ex)
    @eval $ex
end

for fn7 = (:^,), fn8 = (:*,)
    println(fn7)
    ex = quote
        function ($fn7)(a::Dual, b::Float64)
            xx = ($fn7)(a.fun, b)
            yy = ($fn8)(b,($fn7)(a.fun,b-1),a.der)
            return Dual(xx, yy)
        end
    end
    println(ex)
    @eval $ex
end

for fn9 = (:+, :-)
    println(fn9)
    ex = quote
        function ($fn9)(a::Dual)
            xx = ($fn9)(a.fun)
            yy = ($fn9)(a.der)
            return Dual(xx, yy)
        end
    end
    println(ex)
    @eval $ex
end


#Aqui definimos lo que hace falta para la multiplicación y division de un número con un dual, y la suma y resta 
#con numeros reales.
function *(a::Real,b::Dual)
    Dual(a)*b
end

function *(a::Dual,b::Real)
    return b*a
end

#Para la suma y la resta, al derivar una constante esta se hace cero por lo que la segunda entrada no debe cambiar.
function +(a::Real,b::Dual)
    Dual(a)+b
end

function +(a::Dual,b::Real)
    return b+a
end

function -(a::Real,b::Dual)
    Dual(a)-b
end

function -(a::Dual,b::Real)
    return a-Dual(b)
end

#La división
function /(a::Dual,b::Real)
    a/Dual(b)
end

function /(a::Real,b::Dual)
    Dual(a)/b
end
    


end
