"""
lorenz_all_fast(t0,tf,x0,y0,z0,r,sigma,b,p,TOL)

Resuleve el sistema Lorenz parametrizado con r, sigma, b y condiciones iniciales
x0, y0, z0 integrando desde t0 a tf. p es el orden entero de los Taylors utilizados
en el cálculo y TOL es un  parámetro que indica la convergencia de los expoenentes
de Lyapunov. Devuelve

tl, xl, yl, zl, xtl, ytl, ztl, lambda1l, lambda2l, lambda3l, flag

donde tl es una lista de tiempos entre t0 y tf. xl, yl, zl son las listas que
contienen la trayectoria integrada a cada tiempo de la lista tl. xtl, ytl y ztl
son listas con los `Taylors` de cada tiempo de tl. lambda1l, lambda2l, lambda3l
son listas de los tres exponentes calculados a cada tiempo de tl. flag es un
booleano que indica si los exponentes de Lyapunov convergieron en el tiempo de
integración.

El cálculo de los exponentes se hace por medio de un algoritmo diferencial para
las ecuaciones variacionales que resulta ser mucho más rápido que el típico
proceso de Gramm-Schmidt y con menor fluctuaciones.
"""
function lorenz_all_fast(t0::Real, tf::Real, x0::Real, y0::Real, z0::Real, r::Real, sigma::Real = 10.0, b::Real = 8/3, p::Int = 20, TOL::Real = 1.0e-8)
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

		#Arreglo para la integración del jacobiano
		J_int = zeros(3,3)
		#Arreglo auxiliar para integrar las series de Taylor
		n = Array{Real,1}(p)
		for i in range(1,p)
		    n[i] = 1/i
		end


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

      #A continuación se calcula la serie de Taylor hasta orden p
      #usando las relaciones de recurrencia entre las x,y,z y sus derivadas.
      for i in range(1,p)
          #En cada paso se vuelve a calcular la serie de Taylor de dx/dt = f(x,y,z),
          #para dx/dt, dy/dt, dz/dt cada vez a mayor orden.
          dx = sigma*(y-x)
          dy = r*x - y - x*z
          dz = x*y - b*z
          #De ésta se extrae el nuevo coeficiente de x(t), y(t), z(t) que se anexa.
          x = Taylor(push!(x.taylor_vec,dx.taylor_vec[i]/i))
          y = Taylor(push!(y.taylor_vec,dy.taylor_vec[i]/i))
          z = Taylor(push!(z.taylor_vec,dz.taylor_vec[i]/i))
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
      push!(h,(1/2)*(eps(1.0)/abs(x.taylor_vec[p]))^(1/(p-1)))
      push!(h,(1/2)*(eps(1.0)/abs(y.taylor_vec[p]))^(1/(p-1)))
      push!(h,(1/2)*(eps(1.0)/abs(z.taylor_vec[p]))^(1/(p-1)))
      #Luego se toma h como el mínimo de las anteriores.
      hmin = minimum(h)
      #Aquí hay un condicional importante que permite que el último tiempo sea tf.
      h = (t + hmin < tf) ? hmin : tf - t

      #Crea copias de los Taylors para que cuando se modifiquen en el jacobiano no
      #se modifiquen gloablmente:
      xj = Taylor(copy(x.taylor_vec))
      yj = Taylor(copy(y.taylor_vec))
      xj = Taylor(copy(z.taylor_vec))

			#Integra el jacobiano del tiempo t-h a t:
			J = [-Taylor(sigma,p) Taylor(sigma,p) Taylor(0,p); Taylor(r,p)-zj -Taylor(1,p) -xj; yj xj -Taylor(b,p)];
			for i in 1:3
				for j in 1:3
					J[i,j].taylor_vec[2:p+1] = J[i,j].taylor_vec[1:p] .* n
					J[i,j].taylor_vec[1] = 0
				end
			end

      #Ahora sumamos las series de acuerdo al método de Horner:
      sx = x.taylor_vec[p+1]
      sy = y.taylor_vec[p+1]
      sz = z.taylor_vec[p+1]
			sJ = [J[1,1].taylor_vec[p+1] J[1,2].taylor_vec[p+1] J[1,3].taylor_vec[p+1];
						J[2,1].taylor_vec[p+1] J[2,2].taylor_vec[p+1] J[2,3].taylor_vec[p+1];
						J[3,1].taylor_vec[p+1] J[3,2].taylor_vec[p+1] J[3,3].taylor_vec[p+1]]

      for k in range(1,p)
          sx = x.taylor_vec[p+1-k] + h*sx
          sy = y.taylor_vec[p+1-k] + h*sy
          sz = z.taylor_vec[p+1-k] + h*sz
					sJ = [J[1,1].taylor_vec[p+1-k] J[1,2].taylor_vec[p+1-k] J[1,3].taylor_vec[p+1-k];
								J[2,1].taylor_vec[p+1-k] J[2,2].taylor_vec[p+1-k] J[2,3].taylor_vec[p+1-k];
								J[3,1].taylor_vec[p+1-k] J[3,2].taylor_vec[p+1-k] J[3,3].taylor_vec[p+1-k]] + h*sJ
      end

      #x0 es ahora la suma de la serie: x(t0+h), y análogo para y,z.
      x0 = sx
      y0 = sy
      z0 = sz
			#La integral del jacobiano se va sumando intervalo por intervalo
			J_int += (sJ + sJ')/2
			lambda1, lambda2, lambda3 = eigvals(J_int)/(t-t0)

			#Llena las listas de respuesta, si flag es false se llena con el anterior, claramente.
      push!(lambda1l, lambda1)
      push!(lambda2l, lambda2)
      push!(lambda3l, lambda3)

      #La condición de salida del loop. Nótese que se tiene que formar el último Taylor en tf y el último exponente de Lyapunov
      if t == tf
        break
      end

      #Se suma el paso al tiempo.
      t += h

    end

    #Devuelve todas las diversas listas y la flag para ver si el cálculo de los exponentes convergió
    return tl, xl, yl, zl, xtl, ytl, ztl, lambda1l, lambda2l, lambda3l
end
