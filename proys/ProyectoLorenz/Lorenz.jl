__precompile__(true)
"""
## LZ

Módulo que resuelve las ecuaciones de Lorenz clásicas

``dx/dt = σ(y-x),``
``dy/dt = rx - y - xz,``
``dz/dt = xy -bz.``

Donde r, σ, b son parámetros reales positivos.

La solución se hace por medio del método de Taylor con polinomios de Taylor
de grado indicado por el usuario.

Véase las funciones

```julia
lorenz
diflorenz
taylor_lorenz
jacobian
lorenz_all
conteo_fractal
dim_fractal_cajas
```
"""
module LZ
	export lorenz
	export diflorenz
	export taylor_lorenz
	export jacobian
	export lorenz_all
	export conteo_fractal
	export dim_fractal_cajas
using TS

"""
## Integrador de las ecuaciones de Lorenz

		lorenz(t0,tf,x0,y0,z0,r,σ,b,p)

Resuelve las ecuaciones de Lorenz con parámetros `r`, `σ`, `b` y condiciones
iniciales `(x0,y0,z0)` del tiempo inicial `t0` al tiempo final `tf`. Las solución
se hace mediante el método de Taylor con polinomios de Taylor hasta orden `p`.

Devuelve cuatro listas
`t, x, y, z`
con las soluciones.

### Ejemplo:
Para obtener la solución del sistema con parámetros `r, σ, b = 28.0, 10.0, 8/3`,
del tiempo inicial `t0 = 0.0` al tiempo final `tf = 100.0` con condiciones
iniciales `(x0, y0, z0) = (0.1, 0.1, 0.1)` con series de Taylor a orden `p = 20`.
```julia
t, x, y, z = lorenz(0.0, 100.0, 0.1, 0.1, 0.1, 28.0, 10.0, 8/3, 20)
```
"""
function lorenz(t0::Real, tf::Real, x0::Real, y0::Real, z0::Real, r::Real, σ::Real, b::Real, p::Int)
    # Inicializamos la lista de respuestas.
		tl = Array{Real}(0)
		xl = Array{Real}(0)
		yl = Array{Real}(0)
		zl = Array{Real}(0)

		#Incializa el tiempo
		t = t0

    # Ejecutamos mientras que el tiempo inicial sea menor al final.
    while true
				#Llena el arreglo con la respuesta previa
				push!(tl,t)
				push!(xl,x0)
				push!(yl,y0)
				push!(zl,z0)
				#La condcición de paro/salida del loop
				if t == tf
					break
				end

        # Empezamos los Taylors con la condición inicial.
        x = Taylor(x0)
        y = Taylor(y0)
        z = Taylor(z0)

        # A continuación se calcula la serie de Taylor hasta orden p
        # usando las relaciones de recurrencia entre las x,y,z y sus derivadas.
        for i in range(1,p)
           # En cada paso se vuelve a calcular la serie de Taylor de dx/dt = f(x,y,z),
           # para dx/dt, dy/dt, dz/dt cada vez a mayor orden.
            dx = σ*(y-x)
            dy = r*x - y - x*z
            dz = x*y - b*z
           # De ésta se extrae el nuevo coeficiente de x(t), y(t), z(t) que se anexa.
            x = Taylor(push!(x.taylor_vec,dx.taylor_vec[i]/i))
            y = Taylor(push!(y.taylor_vec,dy.taylor_vec[i]/i))
            z = Taylor(push!(z.taylor_vec,dz.taylor_vec[i]/i))
        end

        # Ahora se escoge un paso. Como se recomienda, se toman los dos
        # últimos términos de la serie para calcular h1 y h2, mismo que
				# se hace para x,y,z.
        hx1 = (1/2)*(eps(1.0)/abs(x.taylor_vec[p+1]))^(1/p)
        hx2 = (1/2)*(eps(1.0)/abs(x.taylor_vec[p]))^(1/(p-1))
        hy1 = (1/2)*(eps(1.0)/abs(y.taylor_vec[p+1]))^(1/p)
        hy2 = (1/2)*(eps(1.0)/abs(y.taylor_vec[p]))^(1/(p-1))
        hz1 = (1/2)*(eps(1.0)/abs(z.taylor_vec[p+1]))^(1/p)
        hz2 = (1/2)*(eps(1.0)/abs(z.taylor_vec[p]))^(1/(p-1))
        # Luego se toma h como el mínimo de las anteriores
				hmin = minimum([hx1, hx2, hy1, hy2, hz1, hz2])
				# Si el paso hace que t se pase de tf, se toma tal que t+h=tf
        h = (t + hmin < tf) ? hmin : tf - t

        # Ahora sumamos las series de acuerdo al método de Horner:
        sx = x.taylor_vec[p+1]
        sy = y.taylor_vec[p+1]
        sz = z.taylor_vec[p+1]
        for k in range(1,p)
            sx = x.taylor_vec[p+1-k] + h*sx
            sy = y.taylor_vec[p+1-k] + h*sy
            sz = z.taylor_vec[p+1-k] + h*sz
        end

        # x0 es ahora la suma de la serie: x(t0+h), y análogo para y,z.
        x0 = sx
        y0 = sy
        z0 = sz

				#Actualiza el tiempo
				t += h

    end
    # Se devuelven las listas.
    return tl, xl, yl, zl
end

"""
## Comparador de Trayectorias

		diflorenz(t0,tf,x01,y01,z01,x02,y02,z02,r,σ,b,p)

Resuelve las ecuaciones de Lorenz con parámetros `r`, `σ`, `b` para dos condiciones
iniciales diferentes. La primera es un sistema con `(x01, y01, z01)` y el segundo
sistema es `(x02, y02, z02)`. Los dos sistemas se resuelven a los mismos intervalos
de tiempo dados en la lista `t`.

Devuelve siete listas
`t, x1, y1, z1, x2, y2, z2`
con las soluciones de los dos sistemas y la lista de tiempos.

### Ejemplo:
```julia
t, x1, y1, z1, x2, y2, z2 = diflorenz(0.0, 100.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 + 1.0e-9, 28.0, 10.0, 8/3, 20)
```
Calcula las trayectorias de `t∈[0,100]` para dos sistemas que difieren en condiciones
iniciales por 1.0e-9 en z. Los parámetros son `r, σ, b = 28, 10, 8/3` y el orden de los
polinomios de Taylor es `p=20`.
"""
function diflorenz(t0::Real, tf::Real, x01::Real, y01::Real, z01::Real, x02::Real, y02::Real, z02::Real, r::Real, σ::Real, b::Real, p::Int)
		# Esta función es prácticamente igual a la anterior, lorenz(...) que calcula una solución.
		# La diferencia central es que el paso de tiempo se selcciona como el mínimo del paso de
		# cada solución para tener soluciones en los mismos tiempos.
		# Inicializamos la lista de respuestas.
		tl = Array{Real}(0)
		x1l = Array{Real}(0)
		y1l = Array{Real}(0)
		z1l = Array{Real}(0)
		x2l = Array{Real}(0)
		y2l = Array{Real}(0)
		z2l = Array{Real}(0)

		#Incializa el tiempo
		t = t0

    while true
				#Llena las listas con la respuesta anterior
				push!(tl,t)
				push!(x1l,x01)
				push!(y1l,y01)
				push!(z1l,z01)
				push!(x2l,x02)
				push!(y2l,y02)
				push!(z2l,z02)
				#La condición de salida del loop
				if t == tf
					break
				end

				x1 = Taylor(x01)
        y1 = Taylor(y01)
        z1 = Taylor(z01)

        x2 = Taylor(x02)
        y2 = Taylor(y02)
        z2 = Taylor(z02)

        for i in range(1,p)
            dx1 = σ*(y1-x1)
            dy1 = r*x1-y1-x1*z1
            dz1 = x1*y1 - b*z1
            dx2 = σ*(y2-x2)
            dy2 = r*x2-y2-x2*z2
            dz2 = x2*y2 - b*z2

            x1 = Taylor(push!(x1.taylor_vec,dx1.taylor_vec[i]/i))
            y1 = Taylor(push!(y1.taylor_vec,dy1.taylor_vec[i]/i))
            z1 = Taylor(push!(z1.taylor_vec,dz1.taylor_vec[i]/i))
            x2 = Taylor(push!(x2.taylor_vec,dx2.taylor_vec[i]/i))
            y2 = Taylor(push!(y2.taylor_vec,dy2.taylor_vec[i]/i))
            z2 = Taylor(push!(z2.taylor_vec,dz2.taylor_vec[i]/i))
        end

        hx11 = (1/2)*(eps(1.0)/abs(x1.taylor_vec[p+1]))^(1/p)
        hx12 = (1/2)*(eps(1.0)/abs(x1.taylor_vec[p]))^(1/(p-1))
        hy11 = (1/2)*(eps(1.0)/abs(y1.taylor_vec[p+1]))^(1/p)
        hy12 = (1/2)*(eps(1.0)/abs(y1.taylor_vec[p]))^(1/(p-1))
        hz11 = (1/2)*(eps(1.0)/abs(z1.taylor_vec[p+1]))^(1/p)
        hz12 = (1/2)*(eps(1.0)/abs(z1.taylor_vec[p]))^(1/(p-1))
        hx21 = (1/2)*(eps(1.0)/abs(x2.taylor_vec[p+1]))^(1/p)
        hx22 = (1/2)*(eps(1.0)/abs(x2.taylor_vec[p]))^(1/(p-1))
        hy21 = (1/2)*(eps(1.0)/abs(y2.taylor_vec[p+1]))^(1/p)
        hy22 = (1/2)*(eps(1.0)/abs(y2.taylor_vec[p]))^(1/(p-1))
        hz21 = (1/2)*(eps(1.0)/abs(z2.taylor_vec[p+1]))^(1/p)
        hz22 = (1/2)*(eps(1.0)/abs(z2.taylor_vec[p]))^(1/(p-1))

				# Esta es la diferencia central a la función anterior, se selecciona
				# la misma h considerando las seis distintas. Así los pasos de tiempo
				# son iguales en las dos soluciones.
        hmin = minimum([hx11, hx12, hy11, hy12, hz11, hz12, hx21, hx22, hy21, hy22, hz21, hz22])
        h = (t + hmin < tf) ? hmin : tf - t

				#Análogamente sumamos por el método de Horner
        sx1 = x1.taylor_vec[p+1]
        sy1 = y1.taylor_vec[p+1]
        sz1 = z1.taylor_vec[p+1]
        sx2 = x2.taylor_vec[p+1]
        sy2 = y2.taylor_vec[p+1]
        sz2 = z2.taylor_vec[p+1]
        for k in range(1,p)
            sx1 = x1.taylor_vec[p+1-k] + h*sx1
            sy1 = y1.taylor_vec[p+1-k] + h*sy1
            sz1 = z1.taylor_vec[p+1-k] + h*sz1
            sx2 = x2.taylor_vec[p+1-k] + h*sx2
            sy2 = y2.taylor_vec[p+1-k] + h*sy2
            sz2 = z2.taylor_vec[p+1-k] + h*sz2
        end

        x01 = sx1
        y01 = sy1
        z01 = sz1
        x02 = sx2
        y02 = sy2
        z02 = sz2

				#Actualiza el tiempo
				t += h

    end
    return tl, x1l, y1l, z1l, x2l, y2l, z2l
end

"""
## Integrador de las ecuaciones de Lorenz

		taylor_lorenz(t0,tf,x0,y0,z0,r,σ,b,p)

Resuelve las ecuaciones de Lorenz con parámetros `r, σ, b` por el método de
Taylor, al igual que `lorenz(t0, tf, x0, y0, z0, r, σ, b, p)`, con la diferencia de que
devuelve

`tl, xl, yl, zl`

Donde `tl` es un arreglo con los tiempos entre `t0` y `tf` en pasos que elige el método
y `xl, yl, zl` son arreglos con los `Taylor`'s correspondientes a cada tiempo.

"""
function taylor_lorenz(t0::Real, tf::Real, x0::Real, y0::Real, z0::Real, r::Real, σ::Real, b::Real, p::Int)
    # Inicializamos la lista de respuestas.
    tl = Array{Real}(0)
    xl = Array{Taylor}(0)
    yl = Array{Taylor}(0)
    zl = Array{Taylor}(0)

    # Ejecutamos mientras que el tiempo inicial sea menor al final.
    t = t0
    while true
        # Empezamos los Taylors con la condición inicial.
        x = Taylor(x0)
        y = Taylor(y0)
        z = Taylor(z0)

        # A continuación se calcula la serie de Taylor hasta orden p
        # usando las relaciones de recurrencia entre las x,y,z y sus derivadas.
        for i in range(1,p)
           # En cada paso se vuelve a calcular la serie de Taylor de dx/dt = f(x,y,z),
           # para dx/dt, dy/dt, dz/dt cada vez a mayor orden.
            dx = σ*(y-x)
            dy = r*x - y - x*z
            dz = x*y - b*z
           # De ésta se extrae el nuevo coeficiente de x(t), y(t), z(t) que se anexa.
            x = Taylor(push!(x.taylor_vec,dx.taylor_vec[i]/i))
            y = Taylor(push!(y.taylor_vec,dy.taylor_vec[i]/i))
            z = Taylor(push!(z.taylor_vec,dz.taylor_vec[i]/i))
        end

				#Llena las listas con los Taylors
        push!(tl,t)
        push!(xl,x)
        push!(yl,y)
        push!(zl,z)
				if t == tf
					break
				end

        # Ahora se escoge un paso. Como se recomienda, se toman los dos
        # últimos términos de la serie para calcular h1 y h2, mismo que
        # se hace para x,y,z.
        hx1 = (1/2)*(eps(1.0)/abs(x.taylor_vec[p+1]))^(1/p)
        hx2 = (1/2)*(eps(1.0)/abs(x.taylor_vec[p]))^(1/(p-1))
        hy1 = (1/2)*(eps(1.0)/abs(y.taylor_vec[p+1]))^(1/p)
        hy2 = (1/2)*(eps(1.0)/abs(y.taylor_vec[p]))^(1/(p-1))
        hz1 = (1/2)*(eps(1.0)/abs(z.taylor_vec[p+1]))^(1/p)
        hz2 = (1/2)*(eps(1.0)/abs(z.taylor_vec[p]))^(1/(p-1))
        # Luego se toma h como el mínimo de las anteriores y se suma
        # al tiempo.
        hmin = minimum([hx1, hx2, hy1, hy2, hz1, hz2])
        h = (t + hmin < tf) ? hmin : tf - t
        t += h

        # Ahora sumamos las series de acuerdo al método de Horner:
        sx = x.taylor_vec[p+1]
        sy = y.taylor_vec[p+1]
        sz = z.taylor_vec[p+1]
        for k in range(1,p)
            sx = x.taylor_vec[p+1-k] + h*sx
            sy = y.taylor_vec[p+1-k] + h*sy
            sz = z.taylor_vec[p+1-k] + h*sz
        end

        # x0 es ahora la suma de la serie: x(t0+h), y análogo para y,z.
        x0 = sx
        y0 = sy
        z0 = sz

    end
    # Se devuelven las lista.
    return tl, xl, yl, zl
end

"""
## Jacobiano de las ecuaciones de Lorenz

		jacobian(t, tl, xl, yl, zl, r, σ, b)

Calcula el jacobiano de las ecuaciones de Lorenz para parámetros `r, σ, b` a
un instante `t` sobre una trayectoria integrada. Las entradas `tl, xl, yl, zl`
son arreglos de tiempo y de `Taylor`'s obtenidos con la función `taylor_lorenz`.
Devuelve la matriz jacobiana:

```julia
[[-σ, σ, 0],
 [r-z(t), -1, -x(t)],
 [y(t), x(t), -b]]
```

### Ejemplo
Primero se calcula la solución entre un intervalo de tiempo `t ∈ [0,100]`, por decir
y con ciertos parámetros y condiciones iniciales.
```julia
tl, xl, yl, zl = taylor_lorenz(0.0, 100.0, 0.1, 0.1, 0.1, 28.0, 10.0, 8/3, 20)
```
Ahora se puede calcular el jacobiano sobre cualquier punto de la trayectoria
parametrizado por un instante de tiempo, por ejemplo `t=5.0`
```julia
jacobian(5.0, tl, xl, yl, zl, 28.0, 10.0, 8/3)
```
"""
function jacobian(t::Real,tl::Array{Real,1},xl::Array{Taylor,1},yl::Array{Taylor,1},zl::Array{Taylor,1},r::Real, σ::Real, b::Real)
    #Debemos primero encontrar en qué radio de convergencia se halla el tiempo t
    for i in range(1,length(tl)-1)
        if (tl[i] <= t) & (t <= tl[i+1])
            count = i

            #Extrae el Taylor correspondiente a este radio de convergencia:
            x,y,z = xl[count], yl[count], zl[count]

            #El paso de tiempo corresponde a t-t[i]
            h = t - tl[count]

            #Suma por el método de Horner:
            sx = x.taylor_vec[x.gr+1]
            sy = y.taylor_vec[x.gr+1]
            sz = z.taylor_vec[x.gr+1]
            for k in range(1,x.gr)
                sx = x.taylor_vec[x.gr+1-k] + h*sx
                sy = y.taylor_vec[x.gr+1-k] + h*sy
                sz = z.taylor_vec[x.gr+1-k] + h*sz
            end

            # xf,yf,zf son los flotantes x(t),y(t),z(t)
            xf = sx
            yf = sy
            zf = sz

            #Devuelve la matriz Jacobiana
            return [-σ σ 0; r-zf -1 -xf; yf xf -b]
            break
        end
    end
end

"""
## Integrador de trayectoria y exponentes de Lyapunov para las ecuaciones de Lorenz

		lorenz_all(t0,tf,x0,y0,z0,r,σ,b,p,TOL,gs_count)

Resuleve el sistema Lorenz parametrizado con `r, σ, b` y condiciones iniciales
`x0, y0, z0` integrando desde `t0 a tf`. `p` es el orden entero de los `Taylor`'s utilizados
en el cálculo y `TOL` es un  parámetro que indica la convergencia de los exponentes
de Lyapunov. `gs_count` es cuántos pasos deben suceder para que se vuelve a ortogonalizar por
GS. Devuelve

`tl, xl, yl, zl, xtl, ytl, ztl, lambda1l, lambda2l, lambda3l, flag`

donde `tl` es una lista de tiempos entre `t0` y `tf`. `xl, yl, zl` son las listas que
contienen la trayectoria integrada a cada tiempo de la lista `tl`. `xtl, ytl y ztl`
son listas con los `Taylor`'s de cada tiempo de `tl`. `lambda1l, lambda2l, lambda3l`
son listas de los tres exponentes calculados a cada tiempo de `tl`. `flag` es un
booleano que indica si los exponentes de Lyapunov convergieron en el tiempo de
integración.

El cálculo de los exponentes se hace por medio de las ecuaciones variacionales,
resolviéndolas por el método de Taylor. Se hace otrogonalización de Gramm-Schmidt
en cada paso para evitar que los eigenvectores se "aplanen" en la dirección de
máximo crecimiento.

### Ejemplo
Suponga que se quiere conocer los exponentes de Lyapunov del sistema Lorenz parametrizado por
`r, σ, b = 28, 10, 8/3`. Entonces se integra desde un tiempo `t0=0` a `tf=1000` y con tolerancia
`TOL = 1.0e-4`. Las condiciones iniciales se establecen como `x0, y0, z0 = 0.1, 0.1, 0.1`. Por lo
tanto, al correr
```julia
sol = lorenz_all(0.0, 1000.0, 0.1, 0.1, 0.1, 28.0, 10.0, 8/3, 20, 1.0e-4)
flag, l1, l2, l3 = sol[end][end], sol[end-1][end], sol[end-2], sol[end-3][end]
```
`flag` nos dice si los exponentes convergieron dentro del rango de toleracia `TOL` y las `l1, l2, l3`
son los valores finales de los exponentes de Lyapunov.
"""
function lorenz_all(t0::Real, tf::Real, x0::Real, y0::Real, z0::Real, r::Real, σ::Real = 10.0, b::Real = 8/3, p::Int = 20, TOL::Real = 1.0e-5, gs_count::Int = 0)
    #Incializa la lista de respuestas: tl = time_list, xl = x_list, xtl = x_taylor_list, etc.
    tl = Array{Real}(0)
    xl = Array{Real}(0)
    yl = Array{Real}(0)
    zl = Array{Real}(0)
    #Las listas de Taylors
    xtl = Array{Taylor}(0)
    ytl = Array{Taylor}(0)
    ztl = Array{Taylor}(0)
    #Las listas de los exponentes de Lyapunov
    lambda1l = Array{Real}(0)
    lambda2l = Array{Real}(0)
    lambda3l = Array{Real}(0)

    #Incializa el tiempo
    t = t0

    #Inicialización de tres vectores ortonormales y los exponentes de Lyapunov
    u1 = [1.0 0.0 0.0].'
    u2 = [0.0 1.0 0.0].'
    u3 = [0.0 0.0 1.0].'
    lambda1 = 0.0
    lambda2 = 0.0
    lambda3 = 0.0

		#Contador de ortogonalización GS
		counter = 0

    #La bandera que nos indica si el cálculo de los exponentes debe conitnuar
    flag = true

    #Se ejecuta el loop principal. La salida del loop está en el "if t==tf" de abajo.
    while true
      #Llena el arreglo de escalares con la respuesta previa (o la condición inicial en la primera iteración)
      push!(tl,t)
      push!(xl,x0)
      push!(yl,y0)
      push!(zl,z0)

      #Empezamos los Taylors con la condición inicial.
      x = Taylor(x0)
      y = Taylor(y0)
      z = Taylor(z0)

      #Esta condición revisa si se debe conitnuar el cálculo de los exponentes o ya convergieron
      if (length(lambda1l) > 10) && flag
				#Revisa si el promedio de los últimos diez resulatdos está dentro de TOL del nuevo
				if (abs(mean(lambda1l[end-10:end-1]) -lambda1l[end]) < TOL) && (abs(mean(lambda2l[end-10:end-1]) -lambda2l[end]) < TOL) && (abs(mean(lambda3l[end-10:end-1]) -lambda3l[end]) < TOL)
          flag = false
        end
      end

      if flag
        #Convertimos los vectores a Taylors
        u1x = Taylor(u1[1])
        u1y = Taylor(u1[2])
        u1z = Taylor(u1[3])

        u2x = Taylor(u2[1])
        u2y = Taylor(u2[2])
        u2z = Taylor(u2[3])

        u3x = Taylor(u3[1])
        u3y = Taylor(u3[2])
        u3z = Taylor(u3[3])
      end

      #A continuación se calcula la serie de Taylor hasta orden p
      #usando las relaciones de recurrencia entre las x,y,z y sus derivadas.
      for i in range(1,p)
          #En cada paso se vuelve a calcular la serie de Taylor de dx/dt = f(x,y,z),
          #para dx/dt, dy/dt, dz/dt cada vez a mayor orden.
          dx = σ*(y-x)
          dy = r*x - y - x*z
          dz = x*y - b*z
          #De ésta se extrae el nuevo coeficiente de x(t), y(t), z(t) que se anexa.
          x = Taylor(push!(x.taylor_vec,dx.taylor_vec[i]/i))
          y = Taylor(push!(y.taylor_vec,dy.taylor_vec[i]/i))
          z = Taylor(push!(z.taylor_vec,dz.taylor_vec[i]/i))

          if flag
            #El jacobiano:
            j = [-σ σ 0; r-z -1 -x; y x -b]

            #Las eccuaciones diferenciales para cada componente de cada vector son:
            du1x = (j[1] * u1x + j[4] * u1y + j[7] * u1z)
            du1y = (j[2] * u1x + j[5] * u1y + j[8] * u1z)
            du1z = (j[3] * u1x + j[6] * u1y + j[9] * u1z)

            du2x = (j[1] * u2x + j[4] * u2y + j[7] * u2z)
            du2y = (j[2] * u2x + j[5] * u2y + j[8] * u2z)
            du2z = (j[3] * u2x + j[6] * u2y + j[9] * u2z)

            du3x = (j[1] * u3x + j[4] * u3y + j[7] * u3z)
            du3y = (j[2] * u3x + j[5] * u3y + j[8] * u3z)
            du3z = (j[3] * u3x + j[6] * u3y + j[9] * u3z)

            #Se va generando la serie de Taylor para cada componente de cada vector:
            u1x = Taylor(push!(u1x.taylor_vec,du1x.taylor_vec[i]/i))
            u1y = Taylor(push!(u1y.taylor_vec,du1y.taylor_vec[i]/i))
            u1z = Taylor(push!(u1z.taylor_vec,du1z.taylor_vec[i]/i))

            u2x = Taylor(push!(u2x.taylor_vec,du2x.taylor_vec[i]/i))
            u2y = Taylor(push!(u2y.taylor_vec,du2y.taylor_vec[i]/i))
            u2z = Taylor(push!(u2z.taylor_vec,du2z.taylor_vec[i]/i))

            u3x = Taylor(push!(u3x.taylor_vec,du3x.taylor_vec[i]/i))
            u3y = Taylor(push!(u3y.taylor_vec,du3y.taylor_vec[i]/i))
            u3z = Taylor(push!(u3z.taylor_vec,du3z.taylor_vec[i]/i))
          end
      end

      #Llena las listas con los Taylors
      push!(xtl,x)
      push!(ytl,y)
      push!(ztl,z)

      #Se incializa un arreglo donde se pondrán los diversos pasos
      h = Array{Real}(0)

      push!(h,(1/2)*(eps(1.0)/abs(x.taylor_vec[p+1]))^(1/p))
      push!(h,(1/2)*(eps(1.0)/abs(y.taylor_vec[p+1]))^(1/p))
      push!(h,(1/2)*(eps(1.0)/abs(z.taylor_vec[p+1]))^(1/p))
      if flag
        push!(h, (1/2)*(eps(1.0)/abs(u1x.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u1y.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u1z.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u2x.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u2y.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u2z.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u3x.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u3y.taylor_vec[p+1]))^(1/p))
        push!(h, (1/2)*(eps(1.0)/abs(u3z.taylor_vec[p+1]))^(1/p))
      end
      #Luego se toma h como el mínimo de las anteriores.
      hmin = minimum(h)
      #Aquí hay un condicional importante que permite que el último tiempo sea tf.
      h = (t + hmin < tf) ? hmin : tf - t
      #Se suma el paso al tiempo.
      t += h

      #Ahora sumamos las series de acuerdo al método de Horner:
      sx = x.taylor_vec[p+1]
      sy = y.taylor_vec[p+1]
      sz = z.taylor_vec[p+1]
      if flag
        s1x = u1x.taylor_vec[p+1]
        s1y = u1y.taylor_vec[p+1]
        s1z = u1z.taylor_vec[p+1]
        s2x = u2x.taylor_vec[p+1]
        s2y = u2y.taylor_vec[p+1]
        s2z = u2z.taylor_vec[p+1]
        s3x = u3x.taylor_vec[p+1]
        s3y = u3y.taylor_vec[p+1]
        s3z = u3z.taylor_vec[p+1]
      end

      for k in range(1,p)
          sx = x.taylor_vec[p+1-k] + h*sx
          sy = y.taylor_vec[p+1-k] + h*sy
          sz = z.taylor_vec[p+1-k] + h*sz
          if flag
            s1x = u1x.taylor_vec[p+1-k] + h*s1x
            s1y = u1y.taylor_vec[p+1-k] + h*s1y
            s1z = u1z.taylor_vec[p+1-k] + h*s1z
            s2x = u2x.taylor_vec[p+1-k] + h*s2x
            s2y = u2y.taylor_vec[p+1-k] + h*s2y
            s2z = u2z.taylor_vec[p+1-k] + h*s2z
            s3x = u3x.taylor_vec[p+1-k] + h*s3x
            s3y = u3y.taylor_vec[p+1-k] + h*s3y
            s3z = u3z.taylor_vec[p+1-k] + h*s3z
          end
      end

      #x0 es ahora la suma de la serie: x(t0+h), y análogo para y,z.
      x0 = sx
      y0 = sy
      z0 = sz

      if flag
				if counter == gs_count
	        #Ortogonalizamos por GS
	        v1 = [s1x s1y s1z].'
	        u1 = v1/vecnorm(v1)
	        v2 = [s2x s2y s2z].'
	        v2 -= vecdot(v2, u1)*u1
	        u2 = v2/vecnorm(v2)
	        v3 = [s3x s3y s3z].'
	        v3 -= vecdot(v3, u1)*u1 + vecdot(v3, u2)*u2
	        u3 = v3/vecnorm(v3)
					counter = 0
				else
					v1 = [s1x s1y s1z].'
					u1 = v1/vecnorm(v1)
					v2 = [s2x s2y s2z].'
					u2 = v2/vecnorm(v2)
					v3 = [s3x s3y s3z].'
					u3 = v3/vecnorm(v3)
					counter += 1
				end
        #Calculamos los exponentes de Lyapunov
        lambda1 = (lambda1*(t-h-t0) + log(vecnorm(v1)))/(t-t0)
        lambda2 = (lambda2*(t-h-t0) + log(vecnorm(v2)))/(t-t0)
        lambda3 = (lambda3*(t-h-t0) + log(vecnorm(v3)))/(t-t0)
      end

      #Llena las listas de respuesta, si flag es false se llena con el anterior, claramente.
      push!(lambda1l, lambda1)
      push!(lambda2l, lambda2)
      push!(lambda3l, lambda3)

      #La condición de salida del loop. Nótese que se tiene que formar el último Taylor en tf y el último exponente de Lyapunov
      if t == tf
        break
      end

    end

    #Devuelve todas las diversas listas y la flag para ver si el cálculo de los exponentes convergió
    return tl, xl, yl, zl, xtl, ytl, ztl, lambda1l, lambda2l, lambda3l, flag
end

"""
## Conteo fractal.

		conteo_fractal(x, r)

Toma una lista `x` que contiene tres arreglos de datos que
corresponden a la solución en el espacio fase de un sistema de ecuaciones diferenciales.
La función identifica el tamaño que la solución ocupa y posteriormente divide este espacio
en `r×r×r` cubos. La función cuenta el número de cubos que poseen al menos un punto de la
solución dentro de ellos.
"""
function conteo_fractal(x::Array{Array{Real,1},1},r::Int64)

    #En caso de que demos una partición menor que cero, atrapamos el error.
    if r<1
        error("La partición del espacio debe ser un número entero")
    end

    #Inicializamos un vector donde guardaremos la cota inferior para cada una de
    #nuestras coordenadas.
    cota_inf=zeros(Real,3);
    #Inicializamos un vector donde guardaremos el tamaño de cada caja en cada dimensión.
    tamano=zeros(Real,3);
    #Inicializamos una matriz donde guardaremos los puntos en que hemos dividido
    #los intervalos.
    division=zeros(Real,3,r+1);
    #Inicializamos una matriz de 3x3x3 donde guardaremos el conteo de cajas que contienen
    #al menos un punto de la solución.
    cajas=zeros(Int64,r,r,r);


    #Guardamos las cotas inferiores y los tamaños de intervalo para cada coordenada.
    for i=1:3
        #Definimos al intervalo más pequeño como el mínimo menos un 0.5% de la distancia
        #respecto al máximo.
        cota_inf[i]=minimum(x[i])-0.005*(abs(maximum(x[i])-minimum(x[i])));
        #El tamaño total del intervalo que se estudiará se toma como 101% de la distancia
        #entre mínimo y máximo de los datos.
        tamano[i]=(abs(maximum(x[i])-minimum(x[i])))*(1.01);
    end

    #Guardamos los puntos en que hemos dividido los intervalos.
    for i=1:3
        division[i,1] = cota_inf[i]
        for j = 2:r+1
            division[i,j] = cota_inf[i] + abs((j-1)*tamano[i]/r);
        end
    end

    # Revisamos la caja correspondiente a cada uno de los elementos de la matriz de
    # 3x3x3
    for i = 1:length(x[1])
        for j=2:r+1,k=2:r+1,l=2:r+1
            if  (  (x[1][i]<division[1,j] && x[1][i]>=division[1,j-1])
                && (x[2][i]<division[2,k] && x[2][i]>=division[2,k-1])
                && (x[3][i]<division[3,l] && x[3][i]>=division[3,l-1]))
                if cajas[j-1,k-1,l-1] != 0
                    break;
                else
                    #Si existe al menos un punto de la solución en la caja
                    #, esa caja se "llena"
                    cajas[j-1,k-1,l-1] = 1;
                end
            end
        end
    end

    #Regresamos la cuenta de las cajas que poseen al menos un punto de la curva.
    return sum(cajas)
end

"""
## Dimensión fractal.

		dim_fractal_cajas(x, R)

Toma una lista x que contiene tres arreglos de datos que
corresponden a la solución en el espacio fase de un sistema de ecuaciones diferenciales.
La función utiliza la función `conteo_fractal(x,r)` para dividir el espacio fase
de la solución y contar el número de cajas que contienen al menos un elemento
de la solución en su interior. Este proceso lo realiza para `R` tamaños de caja distintos.
Como salida se obtienen los valores de `r` utilizados, los valores de los conteos para cada
`r` y el valor de la dimensión fractal.
"""
function dim_fractal_cajas{T}(x::Array{Array{T,1},1},R::Int)
    #Se inicializa la salida de la función en ceros.
    salida=zeros(Int64,R);
    for i=1:R
        # Se guardan los valores obtenidos para cada r
        salida[i]=conteo_fractal(x,i)
    end
    #Se regresan los valores
    return R,linreg(log(collect(1:R)),log(salida))[2]
end

# Termina el módulo.
end
