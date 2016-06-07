
__precompile__(true)

using TaylorSeries
using PyPlot
using PyCall
@pyimport matplotlib.animation as anim

function anima3C(c1,c2,c3, nombre::ASCIIString, loop)
    px1 = [x[1] for x in c1[1]]
    px2 = [x[1] for x in c2[1]]
    px3 = [x[1] for x in c3[1]]
    py1 = [x[2] for x in c1[1]]
    py2 = [x[2] for x in c2[1]]
    py3 = [x[2] for x in c3[1]]


    fig = figure(figsize=(10,10))
    cuadros = [[plot(px1[i],py1[i], ",", px2[i],py2[i],  ",", px3[i],py3[i], marker = ".", color ="r")] for i=1:loop:length(px1)]

    animacion = anim.ArtistAnimation(fig, cuadros, interval=200, blit=true)
    animacion[:save](nombre*".mp4", extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end

function muestra_animacion(nombre::ASCIIString)
    display("text/html", string("""<video autoplay controls><source src="data:video/x-m4v;base64,""",base64(open(readbytes,nombre*".mp4")),"""" type="video/mp4"></video>"""))
end


#Agrego Línea para instalar modulo de Taylor (cortesia de Luis)
#Pkg.add("TaylorSeries")

const epsilon = 1.0e-20
const G = 1
function paso_int{T<:Real}(x_0::Taylor1{T})
    orden = x_0.order
    h1 = (epsilon/abs(x_0.coeffs[orden + 1]))^(1/orden)
    h2 = (epsilon/abs(x_0.coeffs[orden]))^(1/(orden - 1))
    min(h1, h2)
end

function Horner{T<:Real, S<:Real}(x_0::Taylor1{S}, h::T)
    n = x_0.order
    suma = zeros(n)
    suma[1] = x_0.coeffs[n]
    for j in 2:n
        suma[j] = x_0.coeffs[n + 1 - j] + h*suma[j - 1]
    end
    suma[n]
end

#function Integrador_1(pos_iniciales, vel_ini, masas::Array{Float64,1}, t0::Float64, tf::Float64, p::Int)
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

    while t <= tf && h > 1e-8
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
        arr_vz1 = Float64[p1[end][3]]
        arr_vz2 = Float64[p2[end][3]]
        arr_vz3 = Float64[p3[end][3]]

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


r12 = ((Taylor_arr_x1 - Taylor_arr_x2)^2 + (Taylor_arr_y1 - Taylor_arr_y2)^2 + (Taylor_arr_z1 - Taylor_arr_z2)^2)^(1/2)
r13 = ((Taylor_arr_x1 - Taylor_arr_x3)^2 + (Taylor_arr_y1 - Taylor_arr_y3)^2 + (Taylor_arr_z1 - Taylor_arr_z3)^2)^(1/2)
r23 = ((Taylor_arr_x2 - Taylor_arr_x3)^2 + (Taylor_arr_y2 - Taylor_arr_y3)^2 + (Taylor_arr_z2 - Taylor_arr_z3)^2)^(1/2)

            #@show r12, r3, r23

            #Para el cuerpo 1
            fx1 = -G*(m2*(Taylor_arr_x1 - Taylor_arr_x2)/(r12)^2 + m3*(Taylor_arr_x1 - Taylor_arr_x3)/(r13)^2)
            fy1 = -G*(m2*(Taylor_arr_y1 - Taylor_arr_y2)/(r12)^2 + m3*(Taylor_arr_y1 - Taylor_arr_y3)/(r13)^2)
            fz1 = -G*(m2*(Taylor_arr_z1 - Taylor_arr_z2)/(r12)^2 + m3*(Taylor_arr_z1 - Taylor_arr_z3)/(r13)^2)
            #Para el cuerpo 2
            fx2 = -G*(m3*(Taylor_arr_x2 - Taylor_arr_x3)/(r23)^2 + m1*(Taylor_arr_x2 - Taylor_arr_x1)/(r12)^2)
            fy2 = -G*(m3*(Taylor_arr_y2 - Taylor_arr_y3)/(r23)^2 + m1*(Taylor_arr_y2 - Taylor_arr_y1)/(r12)^2)
            fz2 = -G*(m3*(Taylor_arr_z2 - Taylor_arr_z3)/(r23)^2 + m1*(Taylor_arr_z2 - Taylor_arr_z1)/(r12)^2)
            #Para el cuerpo 3
            fx3 = -G*(m1*(Taylor_arr_x3 - Taylor_arr_x1)/(r13)^2 + m2*(Taylor_arr_x3 - Taylor_arr_x2)/(r23)^2)
            fy3 = -G*(m1*(Taylor_arr_y3 - Taylor_arr_y1)/(r13)^2 + m2*(Taylor_arr_y3 - Taylor_arr_y2)/(r23)^2)
            fz3 = -G*(m1*(Taylor_arr_z3 - Taylor_arr_z1)/(r13)^2 + m2*(Taylor_arr_z3 - Taylor_arr_z2)/(r23)^2)

            #Incluimos los nuevos coeficientes
            arr_x1 = push!(arr_x1, Taylor_arr_vx1.coeffs[j]/j)
            arr_x2 = push!(arr_x2, Taylor_arr_vx2.coeffs[j]/j)
            arr_x3 = push!(arr_x3, Taylor_arr_vx3.coeffs[j]/j)
            arr_y1 = push!(arr_y1, Taylor_arr_vy1.coeffs[j]/j)
            arr_y2 = push!(arr_y2, Taylor_arr_vy2.coeffs[j]/j)
            arr_y3 = push!(arr_y3, Taylor_arr_vy3.coeffs[j]/j)
            arr_z1 = push!(arr_z1, Taylor_arr_vz1.coeffs[j]/j)
            arr_z2 = push!(arr_z2, Taylor_arr_vz2.coeffs[j]/j)
            arr_z3 = push!(arr_z3, Taylor_arr_vz3.coeffs[j]/j)
            arr_vx1 =push!(arr_vx1, fx1.coeffs[j]/j)
            arr_vx2 =push!(arr_vx2, fx2.coeffs[j]/j)
            arr_vx3 =push!(arr_vx3, fx3.coeffs[j]/j)
            arr_vy1 =push!(arr_vy1, fy1.coeffs[j]/j)
            arr_vy2 =push!(arr_vy2, fy2.coeffs[j]/j)
            arr_vy3 =push!(arr_vy3, fy3.coeffs[j]/j)
            arr_vz1 =push!(arr_vz1, fz1.coeffs[j]/j)
            arr_vz2 =push!(arr_vz2, fz2.coeffs[j]/j)
            arr_vz3 =push!(arr_vz3, fz3.coeffs[j]/j)
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
        h13 = paso_int(Taylor_arr_vz1)
        h14 = paso_int(Taylor_arr_vz2)
        h15 = paso_int(Taylor_arr_vz3)


        #Elegimos el h más pequeño
        h = min(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15)
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
                    rij = :( (($(n_Tarrxj(j)) - $(n_Tarrxj(ℓ)))^2 + ($(n_Tarryj(j)) - $(n_Tarryj(ℓ)))^2 + ($(n_Tarrzj(j)) - $(n_Tarrzj(ℓ)))^2 )^1/2 )
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
        ϵ1 = m1*(0.5*norm(vc1[j])^2 - m2/norm(pc2[j] - pc1[j]) - m3/norm(pc3[j] - pc1[j]))
        ϵ2 = m2*(0.5*norm(vc2[j])^2 - m1/norm(pc1[j] - pc2[j]) - m3/norm(pc3[j] - pc2[j]))
        ϵ3 = m3*(0.5*norm(vc3[j])^2 - m2/norm(pc2[j] - pc3[j]) - m1/norm(pc1[j] - pc3[j]))
        ϵ[j] = ϵ1 + ϵ2 + ϵ3
    end
    ϵ
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
