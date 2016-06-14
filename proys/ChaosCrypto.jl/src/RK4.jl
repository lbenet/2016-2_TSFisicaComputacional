## Para precompilar las cosas
__precompile__(true)

"""
#Modulo RK

Implementación del método de integración de ecuaciones diferenciales Runge-Kutta. 
Solamente para fines ilustrativos, utilizado para probar
las primeras versiones de nuestras integraciones. Deprecado a favor de método con Taylor.
"""


module RK
    export runge4, integrar
    function runge4(f,x,t,dt)
        k1 = f(x,t)
        k2 = f(x+dt*k1/2, t+dt/2)
        k3 = f(x+dt*k2/2,t+dt/2)
        k4 = f(x+dt*k3,t + dt)
        (k1+2k2+2k3+k4)/6
    end

    function integrar(f, x0, t0, t_final, dt)
        # necesito una function f(x,t)
        # que come un vector x, el tiempo t
        # y regresa un vector (el campo vectorial evaluado en x,t)
        
        tiempos = [t0]
        xs = typeof(x0)[x0]

        #x = copy(x0)
        
        x = x0
        
        for t in t0:dt:t_final
            k = runge4(f,x,t,dt)
            x_nueva = x + dt*k
            
            push!(xs, x)
            push!(tiempos, t)
            
            x = x_nueva
            
        end
        
        xs, tiempos
    end
end
