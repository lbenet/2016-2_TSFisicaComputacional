#Modulo de Integración de N Cuerpos
#El siguiente modulo define las funciones necesarias para realizar la integración de las ecuaciones correspondientes al problema
#de tres cuerpos, además de funciones que ayudan a analizar las soluciones obtenidas
#Autores: Daniel (https://github.com/danmarurr) y Fernanda (https://github.com/FernandaPerez)

__precompile__(true)

using TaylorSeries
using PyPlot
using PyCall
@pyimport matplotlib.animation as anim

"""
'anima3C' Recibe la información del integrador 'Integra3' sobre los tres cuerpos para animar las posiciones, recibe también el nombre que desea darse a la animación, 'loop', recibe información sobre el número de saltos que debe dar entre cada posición.
Requiere la previa instalación de 'ffmpeg'
"""
function anima3C(c1,c2,c3, nombre::ASCIIString, loop::Int)
    px1 = [x[1] for x in c1[1]]
    px2 = [x[1] for x in c2[1]]
    px3 = [x[1] for x in c3[1]]
    py1 = [x[2] for x in c1[1]]
    py2 = [x[2] for x in c2[1]]
    py3 = [x[2] for x in c3[1]]


    fig = figure(figsize=(5,5))
    cuadros = [[plot(px1[i],py1[i], ",", px2[i],py2[i],  ",", px3[i],py3[i], marker = ".", color ="r")] for i=1:loop:length(px1)]

    animacion = anim.ArtistAnimation(fig, cuadros, interval=200, blit=true)
    animacion[:save](nombre*".mp4", extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end


"""
'muestra_animacion' contiene el código necesario para mostrar en el notebook un video .mp4 que se halle en el mismo directorio que el notebook donde se trabaja
"""
function muestra_animacion(nombre::ASCIIString)
    display("text/html", string("""<video autoplay controls><source src="data:video/x-m4v;base64,""",base64(open(readbytes,nombre*".mp4")),"""" type="video/mp4"></video>"""))
end


#Agrego Línea para instalar modulo de Taylor (cortesia de Luis)
#Pkg.add("TaylorSeries")

#Definición de constantes
const epsilon = 1.0e-20
const G = 1

"""
Calcula el valor de 'h' para el siguiente paso de integración, es necesario que reciba el valor del Taylor calculado para el parámetro a integrar
"""
function paso_int{T<:Real}(x_0::Taylor1{T})
    orden = x_0.order
    h1 = (epsilon/abs(x_0.coeffs[orden + 1]))^(1/orden)
    h2 = (epsilon/abs(x_0.coeffs[orden]))^(1/(orden - 1))
    min(h1, h2)
end


"""
Calcula el valor del parámetro de x_{n+1} empleando el método de Horner
"""
function Horner{T<:Real, S<:Real}(x_0::Taylor1{S}, h::T)
    n = x_0.order
    suma = zeros(n)
    suma[1] = x_0.coeffs[n]
    for j in 2:n
        suma[j] = x_0.coeffs[n + 1 - j] + h*suma[j - 1]
    end
    suma[n]
end





"""
'Integrador3' Integra las ecuaciones de los tres cuerpos empleando la ley de la gravitación universal hasta el tiempo 'tf' Usando las aproximaciones de Taylor a orden 'p'. Considerando las condiciones iniciales. 'cond_ini' debe ser un arreglo con la información de cada cuerpo organizada de la siguiente forma:

'cond_ini' = [c1, c2, c3]

donde cn = [masa_n, x0_n, y0_n, z0_n, vx0_n, vy0_n, vz0_n]

Regresa cuatro arreglos con la siguiente información:

ts --> tiempos

c1 --> Arreglo con la información de posiciones y velocidades del cuerpo 1

c2 --> Arreglo con la información de posiciones y velocidades del cuerpo 2

c3 --> Arreglo con la información de posiciones y velocidades del cuerpo 3
"""
function Integrador3(cond_ini, tf::Float64, p=28)

    #Para ahorrar al momento de escribir los nombres se indicarán así: pi donde i es el número de cuerpo,
    # y p inidica si es posición o velocidad
    #Creamos los arreglos con los tiempos y las soluciones
    v1 = typeof(zeros(Float64, 3))[[cond_ini[1][5], cond_ini[1][6], cond_ini[1][7]]]
    v2 = typeof(zeros(Float64, 3))[[cond_ini[2][5], cond_ini[2][6], cond_ini[2][7]]]
    v3 = typeof(zeros(Float64, 3))[[cond_ini[3][5], cond_ini[3][6], cond_ini[3][7]]]
    p1 = typeof(zeros(Float64, 3))[[cond_ini[1][2], cond_ini[1][3], cond_ini[1][4]]]
    p2 = typeof(zeros(Float64, 3))[[cond_ini[2][2], cond_ini[2][3], cond_ini[2][4]]]
    p3 = typeof(zeros(Float64, 3))[[cond_ini[3][2], cond_ini[3][3], cond_ini[3][4]]]
    ts = Float64[0.] #arreglo inicial para los tiempos

    m1 = cond_ini[1][1]
    m2 = cond_ini[2][1]
    m3 = cond_ini[3][1]

    t = 0.
    h = 1

    while t <= tf && h > 1e-16
        #Creo arreglos de cada parámetro
        arr_x1 = Float64[p1[end][1]]
        arr_x2 = Float64[p2[end][1]]
        arr_x3 = Float64[p3[end][1]]
        arr_y1 = Float64[p1[end][2]]
        arr_y2 = Float64[p2[end][2]]
        arr_y3 = Float64[p3[end][2]]
        arr_z1 = Float64[p1[end][3]]
        arr_z2 = Float64[p2[end][3]]
        arr_z3 = Float64[p3[end][3]]
        arr_vx1 = Float64[v1[end][1]]
        arr_vx2 = Float64[v2[end][1]]
        arr_vx3 = Float64[v3[end][1]]
        arr_vy1 = Float64[v1[end][2]]
        arr_vy2 = Float64[v2[end][2]]
        arr_vy3 = Float64[v3[end][2]]
        arr_vz1 = Float64[v1[end][3]]
        arr_vz2 = Float64[v2[end][3]]
        arr_vz3 = Float64[v3[end][3]]

        #Creo Taylor's para cada parámetro
        for j in 1:p
            Taylor_arr_x1 = Taylor1(arr_x1)
            Taylor_arr_x2 = Taylor1(arr_x2)
            Taylor_arr_x3 = Taylor1(arr_x3)
            Taylor_arr_y1 = Taylor1(arr_y1)
            Taylor_arr_y2 = Taylor1(arr_y2)
            Taylor_arr_y3 = Taylor1(arr_y3)
            Taylor_arr_z1 = Taylor1(arr_z1)
            Taylor_arr_z2 = Taylor1(arr_z2)
            Taylor_arr_z3 = Taylor1(arr_z3)
            Taylor_arr_vx1 = Taylor1(arr_vx1)
            Taylor_arr_vx2 = Taylor1(arr_vx2)
            Taylor_arr_vx3 = Taylor1(arr_vx3)
            Taylor_arr_vy1 = Taylor1(arr_vy1)
            Taylor_arr_vy2 = Taylor1(arr_vy2)
            Taylor_arr_vy3 = Taylor1(arr_vy3)
            Taylor_arr_vz1 = Taylor1(arr_vz1)
            Taylor_arr_vz2 = Taylor1(arr_vz2)
            Taylor_arr_vz3 = Taylor1(arr_vz3)

            #@show Taylor_arr_x1
            ##Definimos la operación de las 6 ecs de movimiento.


r12 = ((Taylor_arr_x1 - Taylor_arr_x2)^2 + (Taylor_arr_y1 - Taylor_arr_y2)^2 + (Taylor_arr_z1 - Taylor_arr_z2)^2)^(3/2)
r13 = ((Taylor_arr_x1 - Taylor_arr_x3)^2 + (Taylor_arr_y1 - Taylor_arr_y3)^2 + (Taylor_arr_z1 - Taylor_arr_z3)^2)^(3/2)
r23 = ((Taylor_arr_x2 - Taylor_arr_x3)^2 + (Taylor_arr_y2 - Taylor_arr_y3)^2 + (Taylor_arr_z2 - Taylor_arr_z3)^2)^(3/2)

            #@show r12, r3, r23
            
            #Definimos los vectores unitarios

            #Para el cuerpo 1
            fx1 = -G*(m2*(Taylor_arr_x1 - Taylor_arr_x2)/(r12) + m3*(Taylor_arr_x1 - Taylor_arr_x3)/(r13))
            fy1 = -G*(m2*(Taylor_arr_y1 - Taylor_arr_y2)/(r12) + m3*(Taylor_arr_y1 - Taylor_arr_y3)/(r13))
            fz1 = -G*(m2*(Taylor_arr_z1 - Taylor_arr_z2)/(r12) + m3*(Taylor_arr_z1 - Taylor_arr_z3)/(r13))
            #Para el cuerpo 2
            fx2 = -G*(m3*(Taylor_arr_x2 - Taylor_arr_x3)/(r23) + m1*(Taylor_arr_x2 - Taylor_arr_x1)/(r12))
            fy2 = -G*(m3*(Taylor_arr_y2 - Taylor_arr_y3)/(r23) + m1*(Taylor_arr_y2 - Taylor_arr_y1)/(r12))
            fz2 = -G*(m3*(Taylor_arr_z2 - Taylor_arr_z3)/(r23) + m1*(Taylor_arr_z2 - Taylor_arr_z1)/(r12))
            #Para el cuerpo 3
            fx3 = -G*(m1*(Taylor_arr_x3 - Taylor_arr_x1)/(r13) + m2*(Taylor_arr_x3 - Taylor_arr_x2)/(r23))
            fy3 = -G*(m1*(Taylor_arr_y3 - Taylor_arr_y1)/(r13) + m2*(Taylor_arr_y3 - Taylor_arr_y2)/(r23))
            fz3 = -G*(m1*(Taylor_arr_z3 - Taylor_arr_z1)/(r13) + m2*(Taylor_arr_z3 - Taylor_arr_z2)/(r23))

            #Incluimos los nuevos coeficientes
            arr_vx1 =push!(arr_vx1, fx1.coeffs[j]/j)
            arr_vx2 =push!(arr_vx2, fx2.coeffs[j]/j)
            arr_vx3 =push!(arr_vx3, fx3.coeffs[j]/j)
            arr_vy1 =push!(arr_vy1, fy1.coeffs[j]/j)
            arr_vy2 =push!(arr_vy2, fy2.coeffs[j]/j)
            arr_vy3 =push!(arr_vy3, fy3.coeffs[j]/j)
            arr_vz1 =push!(arr_vz1, fz1.coeffs[j]/j)
            arr_vz2 =push!(arr_vz2, fz2.coeffs[j]/j)
            arr_vz3 =push!(arr_vz3, fz3.coeffs[j]/j)
            
            arr_x1 = push!(arr_x1, Taylor_arr_vx1.coeffs[j]/j)
            arr_x2 = push!(arr_x2, Taylor_arr_vx2.coeffs[j]/j)
            arr_x3 = push!(arr_x3, Taylor_arr_vx3.coeffs[j]/j)
            arr_y1 = push!(arr_y1, Taylor_arr_vy1.coeffs[j]/j)
            arr_y2 = push!(arr_y2, Taylor_arr_vy2.coeffs[j]/j)
            arr_y3 = push!(arr_y3, Taylor_arr_vy3.coeffs[j]/j)
            arr_z1 = push!(arr_z1, Taylor_arr_vz1.coeffs[j]/j)
            arr_z2 = push!(arr_z2, Taylor_arr_vz2.coeffs[j]/j)
            arr_z3 = push!(arr_z3, Taylor_arr_vz3.coeffs[j]/j)
            
        end
        #Hacemos Taylor de todos los arreglos finales
        Taylor_arr_x1 = Taylor1(arr_x1)
        Taylor_arr_x2 = Taylor1(arr_x2)
        Taylor_arr_x3 = Taylor1(arr_x3)
        Taylor_arr_y1 = Taylor1(arr_y1)
        Taylor_arr_y2 = Taylor1(arr_y2)
        Taylor_arr_y3 = Taylor1(arr_y3)
        Taylor_arr_z1 = Taylor1(arr_z1)
        Taylor_arr_z2 = Taylor1(arr_z2)
        Taylor_arr_z3 = Taylor1(arr_z3)
        Taylor_arr_vx1 = Taylor1(arr_vx1)
        Taylor_arr_vx2 = Taylor1(arr_vx2)
        Taylor_arr_vx3 = Taylor1(arr_vx3)
        Taylor_arr_vy1 = Taylor1(arr_vy1)
        Taylor_arr_vy2 = Taylor1(arr_vy2)
        Taylor_arr_vy3 = Taylor1(arr_vy3)
        Taylor_arr_vz1 = Taylor1(arr_vz1)
        Taylor_arr_vz2 = Taylor1(arr_vz2)
        Taylor_arr_vz3 = Taylor1(arr_vz3)


        #Calculamos todas las h's posibles
        h1 = paso_int(Taylor_arr_x1)
        h2 = paso_int(Taylor_arr_x2)
        h3 = paso_int(Taylor_arr_x3)
        h4 = paso_int(Taylor_arr_y1)
        h5 = paso_int(Taylor_arr_y2)
        h6 = paso_int(Taylor_arr_z3)
        h7 = paso_int(Taylor_arr_z1)
        h8 = paso_int(Taylor_arr_z2)
        h9 = paso_int(Taylor_arr_y3)
        h10 = paso_int(Taylor_arr_vx1)
        h11 = paso_int(Taylor_arr_vx2)
        h12 = paso_int(Taylor_arr_vx3)
        h13 = paso_int(Taylor_arr_vy1)
        h14 = paso_int(Taylor_arr_vy2)
        h15 = paso_int(Taylor_arr_vy3)
        h16 = paso_int(Taylor_arr_vz1)    
        h17 = paso_int(Taylor_arr_vz2)
        h18 = paso_int(Taylor_arr_vz3)
            

            
        #Elegimos el h más pequeño
        h = min(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18)
        t += h
        #Calculamos el siguiente paso usando Horner
        x1 = evaluate(Taylor_arr_x1, h)
        x2 = evaluate(Taylor_arr_x2, h)
        x3 = evaluate(Taylor_arr_x3, h)
        y1 = evaluate(Taylor_arr_y1, h)
        y2 = evaluate(Taylor_arr_y2, h)
        y3 = evaluate(Taylor_arr_y3, h)
        z1 = evaluate(Taylor_arr_z1, h)
        z2 = evaluate(Taylor_arr_z2, h)
        z3 = evaluate(Taylor_arr_z3, h)
        vx1 = evaluate(Taylor_arr_vx1, h)
        vx2 = evaluate(Taylor_arr_vx2, h)
        vx3 = evaluate(Taylor_arr_vx3, h)
        vy1 = evaluate(Taylor_arr_vy1, h)
        vy2 = evaluate(Taylor_arr_vy2, h)
        vy3 = evaluate(Taylor_arr_vy3, h)
        vz1 = evaluate(Taylor_arr_vz1, h)
        vz2 = evaluate(Taylor_arr_vz2, h)
        vz3 = evaluate(Taylor_arr_vz3, h)

        #Revisamos que no salgan singularidades


        #Creamos vectores para la posicion y la velocidad de cada cuerpo
        p1_temp = Float64[x1, y1, z1]
        p2_temp = Float64[x2, y2, z2]
        p3_temp = Float64[x3, y3, z3]
        v1_temp = Float64[vx1, vy1, vz1]
        v2_temp = Float64[vx2, vy2, vz2]
        v3_temp = Float64[vx3, vy3, vz3]

        #Agregamos nueva información de tiempos, posiciones y velocidades
        p1 = push!(p1, p1_temp)
        p2 = push!(p2, p2_temp)
        p3 = push!(p3, p3_temp)
        v1 = push!(v1, v1_temp)
        v2 = push!(v2, v2_temp)
        v3 = push!(v3, v3_temp)
        ts = push!(ts, t)
    end
    #Creamos arreglos con toda la información de los cuerpos
    cuerpo1 = typeof(p1)[p1, v1]
    cuerpo2 = typeof(p2)[p2, v2]
    cuerpo3 = typeof(p3)[p3, v3]

    ts, cuerpo1, cuerpo2, cuerpo3
end

"""
'IntegradorN' Integra las ecuaciones de los N cuerpos empleando la ley de la gravitación universal hasta el tiempo 'tf' Usando las aproximaciones de Taylor a orden 'p'. Considerando las condiciones iniciales. 'cond_ini' debe ser un arreglo con la información de cada cuerpo organizada de la siguiente forma:

'cond_ini' = [c1, c2, ... , cn]

donde ci = [masa_i, x0_i, y0_i, z0_i, vx0_i, vy0_i, vz0_i] 

Regresa cuatro arreglos con la siguiente información:

ts --> tiempos

ci --> Arreglo con la información de posiciones y velocidades del cuerpo i
"""
function IntegradorN(cond_ini, tf, p=28)
    #Creamos símbolos y expresiones las variables internas
    #Nombres para variables de posiciones, velocidades y masas i-ésimas
    n_vj(j) = symbol(string("v", j))
    n_pj(j) = symbol(string("p", j))
    n_mj(j) = symbol(string("m", j))
    #Nombres para arreglos de variables de posiciones y velocidades  i-ésimas
    n_arrxj(j) = symbol(string("arr_x", j))
    n_arryj(j) = symbol(string("arr_y", j))
    n_arrzj(j) = symbol(string("arr_z", j))
    n_arrvxj(j) = symbol(string("arr_vx", j))
    n_arrvyj(j) = symbol(string("arr_vy", j))
    n_arrvzj(j) = symbol(string("arr_vz", j))
    #Nombres para Taylor's de de variables de posiciones y velocidades  i-ésimas
    n_Tarrxj(j) = symbol(string("Tarr_x", j))
    n_Tarryj(j) = symbol(string("Tarr_y", j))
    n_Tarrzj(j) = symbol(string("Tarr_z", j))
    n_Tarrvxj(j) = symbol(string("Tarr_vx", j))
    n_Tarrvyj(j) = symbol(string("Tarr_vy", j))
    n_Tarrvzj(j) = symbol(string("Tarr_vz", j))
    #Nombre para los Taylor's calculados con ecuaciones de movimiento
    n_fxj(j) = symbol(string("fx",j))
    n_fyj(j) = symbol(string("fy",j))
    n_fzj(j) = symbol(string("fz",j))
    #Nombre para los nuevos parámetros calculados
    n_x(j) = symbol(string("x", j))
    n_y(j) = symbol(string("y", j))
    n_z(j) = symbol(string("z", j))
    n_vx(j) = symbol(string("vx", j))
    n_vy(j) = symbol(string("vy", j))
    n_vz(j) = symbol(string("vz", j))
    #Nombre para los arreglos con toda la información
    n_c(j) = symbol(string("c", j))

    #Nombre para todas las h
    n_h(j) = symbol(string("h",j))

    #Determinamos el número de cuerpos
    const N = length(cond_ini)

    #Creamos los vectores con las N posiciones y velocidades, además de tener las N masas separadas
    for j in 1:N
        #Determinamos valores iniciales
        mj = cond_ini[j][1]
        pjx = cond_ini[j][2]
        pjy = cond_ini[j][3]
        pjz = cond_ini[j][4]
        vjx = cond_ini[j][5]
        vjy = cond_ini[j][6]
        vjz = cond_ini[j][7]

        #Creamos las expresiones
        e_vj = :(typeof(zeros(Float64, 3))[[$vjx,$vjy,$vjz]])
        e_pj = :(typeof(zeros(Float64, 3))[[$pjx,$pjy,$pjz]])
        e_mj = :($mj)
        expr = :($(n_pj(j)) = $(e_pj))
        expr2 = :($(n_vj(j)) = $(e_vj))
        expr3 = :($(n_mj(j)) = $(e_mj))
        eval(expr)
        eval(expr2)
        eval(expr3)
    end

    #Creamos el arreglo con los tiempos
    ts = Float64[0.0]
    t = 0.0
    h = 1
    while t <= tf && h > 1e-8
        #Primero, debemos crear los 6*N arreglos con la información de cada parámetro
        for j in 1:N
            #Creamos expresiones correspondientes y usamos los datos
            expr = :($(n_arrxj(j)) = Float64[$(n_pj(j))[end][1]])
            expr2 = :($(n_arryj(j)) = Float64[$(n_pj(j))[end][2]])
            expr3 = :($(n_arrzj(j)) = Float64[$(n_pj(j))[end][3]])
            expr4 = :($(n_arrvxj(j)) = Float64[$(n_vj(j))[end][1]])
            expr5 = :($(n_arrvyj(j)) = Float64[$(n_vj(j))[end][2]])
            expr6 = :($(n_arrvzj(j)) = Float64[$(n_vj(j))[end][3]])
            eval(expr)
            eval(expr2)
            eval(expr3)
            eval(expr4)
            eval(expr5)
            eval(expr6)
        end

        #Ahora, debemos calcular los Taylor's de los parámetros usando las ecuaciones de movimiento desde orden 1
        #hasta el orden p
        for k in 1:p
            #Primero, debemos crear los Taylor's para cada parámetro

            for j in 1:N
                #Creamos expresiones correspondientes y usamos los datos para crear Taylor's
                expr = :($(n_Tarrxj(j)) =  Taylor1($(n_arrxj(j))))
                expr2 = :($(n_Tarryj(j)) = Taylor1($(n_arryj(j))))
                expr3 = :($(n_Tarrzj(j)) = Taylor1($(n_arrzj(j))))
                expr4 = :($(n_Tarrvxj(j)) = Taylor1($(n_arrvxj(j))))
                expr5 = :($(n_Tarrvyj(j)) = Taylor1($(n_arrvyj(j))))
                expr6 = :($(n_Tarrvzj(j)) = Taylor1($(n_arrvzj(j))))
                eval(expr)
                eval(expr2)
                eval(expr3)
                eval(expr4)
                eval(expr5)
                eval(expr6)
            end

            #Ahora, debemos calcular para cada cuerpo sus 3 parámetros: x, y, z
            for j in 1:N
                #Primero, averguamos los N-1 cuerpos con los cuales el cuerpo j tiene interacción
                indices = Int64[]
                for r in 1:N
                    if j != r
                        push!(indices, r)
                    end
                end
                #Ahora, sumemos cada fuerza
                expx = :(0)
                expy = :(0)
                expz = :(0)
                for ℓ in indices
                    rij = :( (($(n_Tarrxj(j)) - $(n_Tarrxj(ℓ)))^2 + ($(n_Tarryj(j)) - $(n_Tarryj(ℓ)))^2 + ($(n_Tarrzj(j)) - $(n_Tarrzj(ℓ)))^2 )^3/2 )
                    expx = :(($expx) - ($(n_mj(ℓ)))*($(n_Tarrxj(j)) - $(n_Tarrxj(ℓ)))/($rij))
                    expy = :(($expy) - ($(n_mj(ℓ)))*($(n_Tarryj(j)) - $(n_Tarryj(ℓ)))/($rij))
                    expz = :(($expz) - ($(n_mj(ℓ)))*($(n_Tarrzj(j)) - $(n_Tarrzj(ℓ)))/($rij))
                    #@show expx
                end
                #Multipliquemos el factor común
                expx = :(($expx)*(G)*($(n_mj(j))))
                expy = :(($expy)*(G)*($(n_mj(j))))
                expz = :(($expz)*(G)*($(n_mj(j))))

                #Creamos las expresiones finales
                expr = :($(n_fxj(j)) = $expx)
                expr2 = :($(n_fyj(j)) = $expy)
                expr3 = :($(n_fzj(j)) = $expz)
                eval(expr)
                eval(expr2)
                eval(expr3)

                #Ya calculamos el Taylor para cada parámetro, ahora usemos las relaciones de recurrencia
                #Agreguemos a los arreglos de cada parámetro los nuevos valores
                expr = :($(n_arrxj(j)) = push!($(n_arrxj(j)) , $(n_Tarrvxj(j)).coeffs[$k]/$k) )
                expr2 = :($(n_arryj(j)) = push!($(n_arryj(j)) , $(n_Tarrvyj(j)).coeffs[$k]/$k) )
                expr3 = :($(n_arrzj(j)) = push!($(n_arrzj(j)) , $(n_Tarrvzj(j)).coeffs[$k]/$k) )
                expr4 = :($(n_arrvxj(j)) = push!($(n_arrvxj(j)) , $(n_fxj(j)).coeffs[$k]/$k) )
                expr5 = :($(n_arrvyj(j)) = push!($(n_arrvyj(j)) , $(n_fyj(j)).coeffs[$k]/$k) )
                expr6 = :($(n_arrvzj(j)) = push!($(n_arrvzj(j)) , $(n_fzj(j)).coeffs[$k]/$k) )
                eval(expr)
                eval(expr2)
                eval(expr3)
                eval(expr4)
                eval(expr5)
                eval(expr6)
            end
        end

         #Ahora sí, calculemos los Taylor's de los arreglos finales para los 6N parámetros
        for j in 1:N
            #Creamos expresiones correspondientes y usamos los datos para crear Taylor's
            expr = :($(n_Tarrxj(j)) = Taylor1($(n_arrxj(j))) )
            expr2 = :($(n_Tarryj(j)) = Taylor1($(n_arryj(j))) )
            expr3 = :($(n_Tarrzj(j)) = Taylor1($(n_arrzj(j))) )
            expr4 = :($(n_Tarrvxj(j)) = Taylor1($(n_arrvxj(j))))
            expr5 = :($(n_Tarrvyj(j)) = Taylor1($(n_arrvyj(j))))
            expr6 = :($(n_Tarrvzj(j)) = Taylor1($(n_arrvzj(j))))
            eval(expr)
            eval(expr2)
            eval(expr3)
            eval(expr4)
            eval(expr5)
            eval(expr6)
        end

        hs = Float64[]
        #Calculamos todas las h's y las almacenamos en hs
        for j in 1:N
            expr = :($(n_h(1*j)) = paso_int($(n_Tarrxj(j))) )
            expr2 = :($(n_h(2*j)) = paso_int($(n_Tarryj(j))) )
            expr3 = :($(n_h(3*j)) = paso_int($(n_Tarrzj(j))) )
            expr4 = :($(n_h(4*j)) = paso_int($(n_Tarrvxj(j))) )
            expr5 = :($(n_h(5*j)) = paso_int($(n_Tarrvyj(j))) )
            expr6 = :($(n_h(6*j)) = paso_int($(n_Tarrvzj(j))) )
            eval(expr)
            eval(expr2)
            eval(expr3)
            eval(expr4)
            eval(expr5)
            eval(expr6)
            e_hxj = :(push!($hs, $(n_h(1*j))))
            e_hyj = :(push!($hs, $(n_h(2*j))))
            e_hzj = :(push!($hs, $(n_h(3*j))))
            e_hvxj = :(push!($hs, $(n_h(4*j))))
            e_hvyj = :(push!($hs, $(n_h(5*j))))
            e_hvzj = :(push!($hs, $(n_h(6*j))))
            eval(e_hxj)
            eval(e_hyj)
            eval(e_hzj)
            eval(e_hvxj)
            eval(e_hvyj)
            eval(e_hvzj)
        end
        #Elegimos a la h buena como la mínima de todas las que están en hs
        h = minimum(hs)
        t += h
        #Ahora hay que calcular los próximos parámetros usando Horner o evaluate
        for j in 1:N
            expr = :($(n_x(j)) = evaluate($(n_Tarrxj(j)), $h))
            expr2 = :($(n_y(j)) = evaluate($(n_Tarryj(j)), $h))
            expr3 = :($(n_z(j)) = evaluate($(n_Tarrzj(j)), $h))
            expr4 = :($(n_vx(j)) = evaluate($(n_Tarrvxj(j)), $h))
            expr5 = :($(n_vy(j)) = evaluate($(n_Tarrvyj(j)), $h))
            expr6 = :($(n_vz(j)) = evaluate($(n_Tarrvzj(j)), $h))
            eval(expr)
            eval(expr2)
            eval(expr3)
            eval(expr4)
            eval(expr5)
            eval(expr6)
        end
        #Ahora hay que agregar ésta información a los vectores de posición y velocidad
        for j in 1:N
            e_vj = :($(n_vj(j)) = push!($(n_vj(j)), [$(n_vx(j)),$(n_vy(j)),$(n_vz(j))]))
            e_pj = :($(n_pj(j)) = push!($(n_pj(j)), [$(n_x(j)),$(n_y(j)),$(n_z(j))]))
            eval(e_vj)
            eval(e_pj)
        end
        ts = push!(ts, t)
    end

    #Ahora, sólo queda acomodar toda la información de forma adecuada y hacer que el programa regrese esos datos
    for j in 1:N
        expr = :($(n_c(j)) = typeof($(n_pj(j)))[$(n_pj(j)), $(n_vj(j))])
        eval(expr)
    end
    expr = symbol(string("ts"))
    for j in 1:N
        expr =  :($expr, $(n_c(j)))
    end
    eval(expr)
end

function CM_3(c1, c2, c3, masas)
    #Extraemos la información de posiciones y velocidades
    vc1 = c1[2]
    vc2 = c2[2]
    vc3 = c3[2]
    pc1 = c1[1]
    pc2 = c2[1]
    pc3 = c3[1]
    #Obtenemos el valor de las masas
    m1 = masas[1]
    m2 = masas[2]
    m3 = masas[3]
    M = m1 + m2 + m3
    
    CM = (m1*pc1 + m2*pc2 + m3*pc3)/M
    CM
end


"""
'Energia_3' calcula la energía de un sistema de tres partículas a cada paso de tiempo, recibe la información proporcionada por el integrador 'Integrador3' correspondiente a posiciones y velocidades de los tres cuerpos.
Devuelve un arreglo con la energía total del sistema en cada paso de tiempo.
"""
function Energia_3(c1, c2, c3, masas)
    #Extraemos la información de posiciones y velocidades
    vc1 = c1[2]
    vc2 = c2[2]
    vc3 = c3[2]
    pc1 = c1[1]
    pc2 = c2[1]
    pc3 = c3[1]
    #Obtenemos el valor de las masas
    m1 = masas[1]
    m2 = masas[2]
    m3 = masas[3]

    ϵ = zeros(length(vc1))
    for j in 1:length(ϵ)
        K1 = m1*0.5*norm(vc1[j])^2
        K2 = m2*0.5*norm(vc2[j])^2
        K3 = m1*0.5*norm(vc3[j])^2
        U = -G*(m1*m2/norm(pc2[j] - pc1[j]) + m1*m3/norm(pc1[j] - pc3[j]) + m2*m3/norm(pc3[j] - pc2[j]))
        ϵ[j] = K1 + K2 + K3 + U
    end
    ϵ
end

"""
'Angular_3' calcula el momento angular de un sistema de tres partículas a cada paso de tiempo, recibe la información proporcionada por el integrador 'Integrador3' correspondiente a posiciones y velocidades de los tres cuerpos.
Devuelve un arreglo con el momento angular total del sistema en cada paso de tiempo.
"""
function Angular_3(c1, c2, c3, masas)
    #Extraemos la información de posiciones y velocidades
    vc1 = c1[2]
    vc2 = c2[2]
    vc3 = c3[2]
    pc1 = c1[1]
    pc2 = c2[1]
    pc3 = c3[1]
    #Obtenemos el valor de las masas
    m1 = masas[1]
    m2 = masas[2]
    m3 = masas[3]
    #Obtenemos información del L_T a cada paso del tiempo
    L = typeof(vc1[1])[]
    for j in 1:length(vc1)
        L1 = m1*(pc1[j])×(vc1[j])
        L2 = m2*(pc2[j])×(vc2[j])
        L3 = m3*(pc3[j])×(vc3[j])
        LT = L1 + L2 + L3
        push!(L, LT)
    end
    L
end

"""
'CalculaCM(C1, C2, m3)' obtiene, a partir de la información de dos cuerpos 'C1' y 'C2' la información sobre posiciones y velocidades de un tercer cuerpo, de tal forma que el centro de masa de éste sistema se encuentre en el origen y se quede estático (por lo menos a tiempo inicial) 'C1' y 'C2' deben contener información de la siguiente forma:
C1 = [m1, x1, y1, z1, vx1, vy1, vz1] donde m1 es el valor de la masa; x1, y1, z1 las coordendas de la posición y vx1, vy1, vz1 los valores de la velocidad en componentes cartesianas.

La función regresa un arreglo que puede ser introducido como argumento de la función 'Integrador3'.
"""
function CalculaCM(C1, C2, m3)
    #Identificamos la información de CS:
    #Primero las posiciones
    p1 = Float64[C1[2], C1[3], C1[4]]
    p2 = Float64[C2[2], C2[3], C2[4]]
    p3 = zeros(p2)
    #Ahora las velocidades
    v1 = Float64[C1[5], C1[6], C1[7]]
    v2 = Float64[C2[5], C2[6], C2[7]]
    v3 = zeros(v1)
    #Por último las las masas
    m1 = C1[1]
    m2 = C2[1]
    
    #Por comodidad queremos que el centro de masa esté en el origen, en ese caso:
    p3 = -(m1*p1 + m2*p2)/m3
    
    #Queremos que el centro de masa se quede estático, entonces:
    v3 = -(m1*v1 + m2*v2)/m3
    
    #Ahora, regresemos la información de tal forma que el integrador lo acepte
    C3 = Float64[m3, p3[1], p3[2], p3[3], v3[1], v3[2], v3[3]]
    CS = typeof(C3)[C1, C2, C3]
    
    CS    
end


function Integrador_Restringido(masas, p0, v0, tf::Float64, p=28)



    #Para ahorrar al momento de escribir los nombres se indicarán así: pi donde i es el número de cuerpo,
    # y p inidica si es posición o velocidad
    #Creamos los arreglos con los tiempos y las soluciones
    v1 = typeof(zeros(Float64, 3))[p0]

    p1 = typeof(zeros(Float64, 3))[v0]

    ts = Float64[0.] #arreglo inicial para los tiempos

    m1 = masas[1]
    m2 = masas[2]

    M = m1 + m2
    mu = m2 / M
    alpha = m1 / M

    t = 0.
    h = 1

    while t <= tf && h > 1e-8
        #Creo arreglos de cada parámetro
        arr_x1 = Float64[p1[end][1]]
        arr_y1 = Float64[p1[end][2]]
        arr_z1 = Float64[p1[end][3]]

        arr_vx1 = Float64[v1[end][1]]
        arr_vy1 = Float64[v1[end][2]]
        arr_vz1 = Float64[p1[end][3]]


        #Creo Taylor's para cada parámetro
        for j in 1:p
            Taylor_arr_x1 = Taylor1(arr_x1)
            Taylor_arr_y1 = Taylor1(arr_y1)
            Taylor_arr_z1 = Taylor1(arr_z1)

            Taylor_arr_vx1 = Taylor1(arr_vx1)
            Taylor_arr_vy1 = Taylor1(arr_vy1)
            Taylor_arr_vz1 = Taylor1(arr_vz1)

            #@show Taylor_arr_x1
            ##Definimos la operación de las 6 ecs de movimiento.



            r1 = ((Taylor_arr_x1 - mu)^2 + (Taylor_arr_y1)^2 + (Taylor_arr_vz1)^2)
            r2 = ((Taylor_arr_x1 + alpha)^2 + (Taylor_arr_y1)^2 + (Taylor_arr_vz1)^2)
            #@show r12, r3, r23

            #Para el cuerpo
            fx1 = 2*Taylor_arr_vx1 - (alpha / (Taylor_arr_x1 - mu) - (mu*(Taylor_arr_x1 + alpha)/r1^3))
            fy1 = -2*Taylor_arr_vy1 + (1 - (alpha / r1^3))*Taylor_arr_y1
            fz1 = - ((alpha/(r1^3)) + (mu / (r2^3)))*Taylor_arr_z1

            #Incluimos los nuevos coeficientes
            arr_x1 = push!(arr_x1, Taylor_arr_vx1.coeffs[j]/j)
            arr_y1 = push!(arr_y1, Taylor_arr_vy1.coeffs[j]/j)
            arr_z1 = push!(arr_z1, Taylor_arr_vz1.coeffs[j]/j)

            arr_vx1 =push!(arr_vx1, fx1.coeffs[j]/j)
            arr_vy1 =push!(arr_vy1, fy1.coeffs[j]/j)
            arr_vz1 =push!(arr_vz1, fz1.coeffs[j]/j)

        end
        #Hacemos Taylor de todos los arreglos finales
        Taylor_arr_x1 = Taylor1(arr_x1)
        Taylor_arr_y1 = Taylor1(arr_y1)
        Taylor_arr_z1 = Taylor1(arr_z1)

        Taylor_arr_vx1 = Taylor1(arr_vx1)
        Taylor_arr_vy1 = Taylor1(arr_vy1)
        Taylor_arr_vz1 = Taylor1(arr_vz1)



        #Calculamos todas las h's posibles
        h1 = paso_int(Taylor_arr_x1)

        h2 = paso_int(Taylor_arr_y1)

        h3 = paso_int(Taylor_arr_z1)
        h4 = paso_int(Taylor_arr_z2)

        h5 = paso_int(Taylor_arr_vx1)

        h6 = paso_int(Taylor_arr_vz1)



        #Elegimos el h más pequeño
        h = min(h1, h2, h3, h4, h5, h6)
        t += h
        #Calculamos el siguiente paso usando Horner
        x1 = evaluate(Taylor_arr_x1, h)

        y1 = evaluate(Taylor_arr_y1, h)

        z1 = evaluate(Taylor_arr_z1, h)

        vx1 = evaluate(Taylor_arr_vx1, h)

        vy1 = evaluate(Taylor_arr_vy1, h)

        vz1 = evaluate(Taylor_arr_vz1, h)


        #Revisamos que no salgan singularidades


        #Creamos vectores para la posicion y la velocidad de cada cuerpo
        p1_temp = Float64[x1, y1, z1]

        v1_temp = Float64[vx1, vy1, vz1]


        #Agregamos nueva información de tiempos, posiciones y velocidades
        p1 = push!(p1, p1_temp)

        v1 = push!(v1, v1_temp)

        ts = push!(ts, t)
    end
    #Creamos arreglos con toda la información de los cuerpos
    cuerpo1 = typeof(p1)[p1, v1]


    ts, cuerpo1
end
