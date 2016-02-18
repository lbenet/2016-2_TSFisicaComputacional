# Conflictos
## ¿Qué pasa?
1. Usuario A y Usuario B cuentan con un archivo inicial `frases.md`.
2. Usuario A agrega una frase a `frases.md` y hace un `git push`.
3. Usuario B, sin conocer los cambios del usuario A, agrega otra frase
   a `frases.md`.
4. Usuario B hace un `git commit` para tomar en cuenta su edición.
5. Usuario B hace un `git pull` y observa un error que le dice que su
   versión de `frases.md` tiene un conflicto con el que está
   _jalando_.
## ¿Cómo solucionarlo?
1. Usuario B puede abrir el archivo y encontrar que en el aparecen,
   mezcladas en cierto formato, las dos versiones.
2. Usuario B modifica el archivo tomando en cuenta los cambios de usuario A y los
   suyos.
3. Hace un `git commit` para tomar en cuenta los cambios.
