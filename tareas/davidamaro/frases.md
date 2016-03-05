We now see that we have to change the principle of simplicity into a
principle of mathematical beauty.

Este libro está dedicado,
con respeto y admiración,
al principio de mínima acción.

## Conflictos

Cuando existió el conflicto de archivos cuando abro el archivo con
`vim` y creo la versión que tenga los dos archivos.

```
$ vim frases.md
$ git add *
$ git commit -m "Se arreglaron los conflictos."
$ git push origin master
```
# Errores

Cuando se tienen conflictos en un archivo, al abrirlo con un editor de
texto el contenido muestra las dos _versiones_ con cierto formato. Al
modificar (y guardarlo) el archivo de la forma en la que deseamos que sean
__mezclados__ podemos hacer un nuevo _commit_ para tomar en cuenta
estos cambios y que esta sea la versión final.
