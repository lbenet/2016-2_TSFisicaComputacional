#22/03/16

#Este módulo implementa la definición del tipo dual, así como la de sus operaciones.

#Rodríguez Rosenblueth César Daniel
#Santiago Santos Eva Yazmín
#Hernández de la Vega Alejandro

__precompile__(true)

module AD
    import Base: +, -, *, /, ^, ==

    export Dual, xdual

#Se define el Dual con sus dos entradas, la función y la derivada
type Dual{D <: Real}
    fun::D
    der::D
    
end

#Definimos el Dual de una constante y promovemos las entradas.
Dual(fun,der) = Dual(promote(fun, der)...)

Dual(a::Real) = Dual(a,0...)

#Definimos una función para identificar la variable independiente. 
function xdual(x)
    if typeof(x) <: Real
        fun=x
        der=1
        return Dual(fun,der)
    end   
end

import Base: +, -, *, /, ^, ==

# Aqui se implementan las operaciones entre Duales y con Duales y reales.

for fn1 = (:+, :-)
    ex = quote
        function ($fn1)(a::Dual, b::Dual)
            xx = ($fn1)(a.fun, b.fun)
            yy = ($fn1)(a.der, b.der)
            return Dual(xx, yy)
        end
    end
    @eval $ex
end

for fn2 in (:*,), fn3 in (:+,)
    ex = quote
        function ($fn2)(a::Dual, b::Dual)
            xx = ($fn2)(a.fun, b.fun)
            yy = ($fn3)(($fn2)(a.fun, b.der),($fn2)(a.der, b.fun))
            return Dual(xx, yy)
        end
    end
    @eval $ex
end


for fn4 = (:/,), fn5 = (:*,), fn6 = (:-,)
    ex = quote
        function ($fn4)(a::Dual, b::Dual)
            xx = ($fn4)(a.fun, b.fun)
            yy = ($fn4)(($fn6)(a.der, ($fn5)(($fn4)(a.fun,b.fun),b.der)),b.fun)
            return Dual(xx, yy)
        end
    end
    @eval $ex
end

for fn7 = (:^,), fn8 = (:*,)
    ex = quote
        function ($fn7)(a::Dual, b::Float64)
            xx = ($fn7)(a.fun, b)
            yy = ($fn8)(b,($fn7)(a.fun,b-1),a.der)
            return Dual(xx, yy)
        end
    end
    @eval $ex
end

for fn9 = (:+, :-)
    ex = quote
        function ($fn9)(a::Dual)
            xx = ($fn9)(a.fun)
            yy = ($fn9)(a.der)
            return Dual(xx, yy)
        end
    end
    @eval $ex
end

#Finalmente definimos las operaciones de Duales con reales.
for fn10 = (:+, :-, :*, :/,)
    ex = quote
        function ($fn10)(a::Real, b::Dual)
            return ($fn10)(Dual(a),b)
        end
        function ($fn10)(a::Dual, b::Real)
            return ($fn10)(a,Dual(b))
        end
    end
    @eval $ex
end
 
#Identificación
        ==(a::Dual, b::Dual) = (a.fun == b.fun && a.der == b.der) ? true : false 


end
