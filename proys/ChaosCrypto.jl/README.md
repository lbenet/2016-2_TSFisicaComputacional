# ChaosCrypto.jl

Uso del caos para transmitir mensajes encriptados

#### Autores: 
- [Yuriko Yamamoto](https://github.com/Yuriyama "Yuriyama")
- [David-Amaro Alcalá](https://github.com/davidamaro "davidamaro")
- [Ignacio Vargas](https://github.com/ignacio-vc "ignacio-vc")

#### Requisitos:
Incluir la siguiente linea al archivo `.juliarc.jl`, localizado en `/home/username/`
```
push!(LOAD_PATH, "/home/username/path/to/module") 
```

Tal que el directorio clonado localmente, `ChaosCrypto.jl`, se encuentre ahi.

Por ejemplo
```
push!(LOAD_PATH, "/home/ignacio/Documents/Physics/2016-2") 
```
con el directorio `ChaosCrypto.jl` dentro de `2016-2`.

Esto debe hacerse antes de correr Julia, o reiniciar Julia despues de realizar el cambio. Si no existe el archivo,
nada mas se crea.

#### Estructura de Directorios:
- docs: Documentacion y Presentacion
- src: modulos .jl, esta bien hacer merge de estos
- test: Tests
- wip: Work In Progress, Jupyter Notebooks y demas

#### Referencias:

##### Caos, Criptologia

[Nonlinear Dynamics and Chaos - Strogatz (2014)](http://libgen.io/get.php?md5=93608D1E7D48FF61D25173674AF85BD7&key=ALGOHY9BRV5DVM5D "Lib Genesis")

[Synchronization of Lorenz-Based Chaotic Circuits with Applications to Communications](http://www.rle.mit.edu/dspg/documents/SynchroofLorenz.pdf "Artículo")

[Chaos Applications in Telecommunications](http://libgen.io/get.php?md5=0C28EB7B594F94B10BDD9B9391228D85&key=OCJBN9OTSRT306XH "Lib Genesis")

[Handbook of Chaos Control](http://libgen.io/get.php?md5=97455994EC81072A20A21293532926D1&key=I33OJ1BYKTFE0R1S "Lib Genesis")

[Using Chaos to Send Secret Messages](http://bulldog2.redlands.edu/fac/joanna_bieri/nonlinear/Chotic_Messages.pdf "powerpoint")

[Edward Lorenz's Strange Attraction](https://logicaltightrope.com/2013/08/29/edward-lorenzs-strange-attraction/ "blog")

[Sending Your Secrets Safely with Chaos](https://logicaltightrope.com/2013/09/01/sending-your-secrets-safely-with-chaos/ "blog")

[A secret message from another dimension](https://web.archive.org/web/20150214122103/http://jellymatter.com/2012/01/04/a-secret-message-from-another-dimension/ "ejemplo")

[The Chaos Cookbook, A Practical Programming Guide](http://libgen.io/get.php?md5=20F14D04E0992220B0093F9F5D3A7551&key=H8DQQ2F1HMPVXGH8 "Lib Genesis")

##### Julia + Audio
[A really brief introduction to audio signal processing in Julia](http://www.seaandsailor.com/audiosp_julia.html "Audio")

##### Cellular Automata (Imagenes)

[Chaotic Encryption Method Based on Life-like Cellular Automata](http://arxiv.org/pdf/1112.6326v1.pdf "articulo")

[The Research of Image Encryption Algorithms Based on Chaos Cellular Automata](https://pdfs.semanticscholar.org/bff7/e1fc9a4201e9b50b16314ceffd13c024edf4.pdf "articulo")

[New Possiblities for Cellular Automata in Cryptography](http://www.criptored.upm.es/cibsi/cibsi2011/info/Ponencias/5.%20New%20Possibilities%20for%20Cellular%20Automata%20in%20Cryptography.pdf "presentacion")


