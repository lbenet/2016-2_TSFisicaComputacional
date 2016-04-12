# Proyectos de fin de semestre

**Fecha importantes:**

Viernes 3 de junio: Envío de proyectos y presentaciones.

Martes 7 de junio: Examen final (presentación de proyectos), a la 1 PM, en el Aula de Cómputo III.


=====================

## Reglas

- Los proyectos **tienen** que hacerse  en equipo, mínimo dos personas, máximo tres.

- Cada equipo debe escoger un tema (ver la lista que aparece abajo), de tal manera que **no** haya dos equipos que desarrollen el mismo tema. La idea es desarrollar algún aspecto particular del tema que sea relevante en el contexto de sistemas dinámicos. Los detalles se discutirán *en clase*, semana a semana.

- Los equipos pueden proponer otros temas distintos. En este caso, deben presentar la idea de lo que quieren hacer de manera escrita, en una cuartilla como máximo, y el proyecto debe ser aceptado. El tema debe estar relacionado con sistemas dinámicos y deberá mostrar la aplicación los métodos vistos en clase.

- Durante la elaboración del proyecto la asistencia a clase *es importante*: cada semana, iremos monitoreando los avances que se van haciendo.

- Los proyectos serán albergados en una cuenta de GitHub de alguno de los miembros del equipo. Es obligación comunicar el repositorio donde se encuentra el proyecto; fecha límite: martes 26 de abril.

- Al finalizar se deberá enviar un ipynb con todos los detalles del proyecto. En base a éste se puede basar la presentación oral, que a priori puede ser en un archivo independiente. El mismo sistema usado para el envío de tareas se usará para los proyectos, en este caso usando el subdirectorio "proys". Los proyectos deben estar documentados adecuadamente.

=====================

## Temas del proyecto

Los temas que se sugieren a continuación involucran la integración de sistemas de ecuaciones diferenciales ordinarias de diversa dimensionalidad. Algunas preguntas interesantes son si el sistema exhibe o no caos, si hay cambios de estabilidad de las soluciones más sencillas (bifurcaciones de puntos fijos, órbitas periódicas, ciclos límite, etc). Los aspectos que se estudien deben ilustrarse a través de cálculos numéricos. Una buena manera de obtener resultados respecto a los exponentes de Liapunov es integrando las llamadas ecuaciones variacionales. En otras ocasiones conviene ilustrar la estructura del espacio fase a través de secciones de Poincaré o mapeos estroboscópicos.

1. Integración del problema gravitacional de 3 (o más cuerpos). En particular, cálculo de órbitas periódicas u órbitas especiales (como por ejemplo las *coreografías*) y su estabilidad. Ver [1].

- Problema restringido de tres cuerpos: estudiar la transferencia de un satélite en una órbita inicialmente cercana a la Tierra a una alrededor de la Luna. Ver [2].

- Dinámica de las ecuaciones de Lorentz clásicas: exponentes de Liapunov y exploración de los parámetros. Ver [3], problemas 9.3.2-7, 9.3.9-10, 9.5.1-3.

- Uso del caos para transmitir mensajes encriptados. Ver [3], pp 335, y los problemas 9.6.3-5.

- Resonancias en un sistema Hamiltoniano periódicamente forzado (osciladores acoplados forzados). Ver [4].

- Estudiar el modelo Kermack-McKendrick de evolución de una epidemia en términos de sus parámetros, o alternativas (ecuaciones de reacción-difusión). En particular, caracterizar los puntos fijos, su estabilidad. Ver [5], capítulo 19.

- Integración de ecuaciones diferenciales con demora, e.g. $\dot{x}(t) = sin x(t-τ)$, con $τ\geq 0$ constante. Ver [6].

- Exploración numérica de la resonacia estocástica clásica: partícula en un pozo doble forzado con ruido. Ver [7].

- Estudiar el movimiento de varias partículas en un círculo, con interacciones de largo alcance. Ver [8].

- Oscilaciones en sistemas biológicos: por ejemplo, modelos de dinámica poblacional, epidemiología o inmunología. Ver [5]?

=====================

## Presentación del proyecto

La presentación del proyecto constituye el examen final del curso (70% de la calificación), y consiste en un seminario de *máximo* 25 minutos y 5 para preguntas. Es necesario tener aceptadas 75% de las tareas.

=====================

## Referencias

[1] C. Simó, New Families of Solutions in N-Body Problems, in C. Casacuberta et al (Eds), European Congress of Mathematics Vol 201, Series Progress in Mathematics (2001), pp 101-115.
http://dx.doi.org/10.1007/978-3-0348-8268-2_6.

[2] Wang Sang Koon, Martin W. Lo, Jerrold E. Marsden and Shane D. Ross, Heteroclinic connections between periodic orbits and resonance transitions in celestial mechanics, Chaos 10, (2000), 427. http://dx.doi.org/10.1063/1.166509.

[3] Steven H. Strogatz, Nonlinear dynamics and chaos (1994).

[4] C. Simó, Global dynamics and fast indicators, in Broer et al (Eds), Global Analysis of Dynamical Systems (2001), pp 373-389.

[5] J. Murray, Mathematical biology (1989).

[6] J.C. Sprott, A simple chaotic delay differential equation, Phys. Letts. A 366 (2007), 397–402.
http://sprott.physics.wisc.edu/pubs/paper304.pdf.

[7] Th. Wellens, V. Shatokhin and Andreas Buchleitner,
Stochastic resonance, Rep. Prog. Phys. 67 (2004), 45.
http://stacks.iop.org/RoPP/67/45

[8] M. Antoni, S. Ruffo, Clustering and relaxation in hamiltonian long-range dynamics, Phys. Rev. E 52 (1995), 2361–2374.
V. Latora, A. Rapisarda, S. Ruffo,
Lyapunov instability and finite size effects in a system with long-range forces, Phys. Rev. Lett. 80 (1998) 692-695. http://arxiv.org/abs/1210.6316
