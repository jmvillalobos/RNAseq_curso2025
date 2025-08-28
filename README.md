# RNAseq_curso2025
Este es el manual para el curso de genómica funcional de 2025
<img width="700" height="700" alt="output" src="https://github.com/user-attachments/assets/bb1e9b8d-e206-4612-87bf-5d57659ed334" />




## Análisis de expresión diferencial
Ahora veremos cómo se realiza un análisis de expresión diferencial. Para ello, utilizaremos la tabla de datos de un experimento realizado por el Dr. Sergio Campos de LANGEBIO. En este trabajo, se pretendía evaluar la respuesta transcripcional de la levadura Saccharomyces cerevisiae ante la carencia de nitrógeno, y además se deseaba investigar si el factor de transcripción Ste12 tiene un papel relevante en dicha respuesta transcripcional. Con este fin, se ideó un experimento de RNA-seq tanto en un medio rico como en uno carente de cualquier fuente de nitrógeno, utilizando tanto la cepa silvestre de S. cerevisiae como la mutante nula para el factor de transcripción.


## Comenzaremos viendo como instalar los programas que necesitamos
```r

#instalando paqueterias de bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")

#necesitamos cargar edgeR y limma, para eso usamos la funsión library
library(limma)
library(edgeR)

```
    Bioconductor version 3.19 (BiocManager 1.30.23), R 4.4.1 (2024-06-14)
    Warning: package(s) not installed when version(s) same as or greater than current; use
    `force = TRUE` to re-install: 'edgeR'
    Old packages: 'nlme'


    
Para comenzar a trabajar con nuestros datos de expresión debemos posicionarnos en el directorio donde esta nuestra tabla de datos, para posteriormente cargar la tabla en R.

Llamaremos a nuestra tabla y comenzaremos a manipularla


```r

# así podemos saber en que carpeta estamos:
getwd()

# en R tambien podemos crear carpetas
dir.create("RNAseq")


# Podemos saber que archivos están en esta dirección usando el siguiente comando: 
list.files()

# nos vemos a otro directorio:
setwd("/cloud/project/RNAseq")

```



## Para crear un archivo donde se almacenarán los resultados del análisis:

```r

outpath = "/cloud/project/RNAseq/results/"

dir.create(outpath, showWarnings=FALSE)

```


## Aquí llamaremos a nuestra tabla de datos y la guardaremos en una variable dentro de R


### Podemos saber cómo luce nuestra tabla cargada y que dimenciones tiene este objeto. 

```r

counts <- read.table("https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/main/Saccharomyces.txt",
                     header = TRUE, row.names = 1, sep = "\t")
# Podemos saber cómo luce nuestra tabla cargada y qué dimensiones tiene este objeto. 
head(counts)

dim(counts)
```
Algo que regularmente queremos hacer al iniciar un análisis de expresión génica, es filtrar los datos por la cantidad mínima de lecturas asociadas a cada gen, a fin de no tener muchos genes con muy baja expresión y que estos nos causen ruido en el análisis.

A manera de ejemplo, aquí les presento un filtrado usando las siguientes instrucciones


```r
# genes que tienen menos de 2 reads y en la suma de todas las filas da menos de 10 lecturas se elimnarán. 

counts = counts[rowSums(cpm(counts) >= 2) >=10,]

# Aquí vemos el uso del comando cpm y nos ayuda a tener los datos de expresión en cuentas por millón.

# despues de aplicar este filtro, nuestra tabla ha cambiado y podemos volver a evaluarla
head(counts)
colnames(counts)
dim(counts)

```



## Además teniendo los datos precargados, podemos visualizarlos en alguna gráfica como lo hemos visto antes
```r
plot(log2(counts[,c("wt_sc_1", "st_sc_1")]), col="gray") #aquí estamos comparando únicamente dos librerías, ustedes puedes probar otras.
```
<img width="710" height="399" alt="Captura de pantalla 2025-08-27 a la(s) 19 01 05" src="https://github.com/user-attachments/assets/392edd73-60e8-4d61-8815-675a708415e6" />


## Ejercicio
¿Recuerdan cómo pedir específicamente algunas columnas? Ahora comenzaremos a usar funciones de edgeR, este programa nos ayudará a determinar los genes diferencialmente expresados.

*Si desean conocer más acerca de este programa pueden consultar el manual desde R, usando la siguiente instrucción.

edgeRUsersGuide()

## Debemos darle formato a nuestra tabla para que edgeR pueda funcionar sobre ella.

```r

# esto es para declarar los grupos que usaremos en el objeto DGEList (objeto de edgeR)
grp = sub("..$", "", colnames(counts)) 


# vean el objeto
grp
# haciendo el objeto DGEList
dge = DGEList(counts=counts, group=grp)
dge
```

<img width="934" height="476" alt="Captura de pantalla 2025-08-27 a la(s) 19 04 35" src="https://github.com/user-attachments/assets/8b6c549d-3bc1-45b4-ada4-97af930b599c" />

```r
#Podemos visualizar la dispersión de nuestras librerías usando el plot MDS

plotMDS(dge)

```
<img width="739" height="448" alt="Captura de pantalla 2025-08-27 a la(s) 19 06 06" src="https://github.com/user-attachments/assets/30fc48fc-5fdb-4605-a1a3-8bef2d5d2acc" />

## ¿Réplicas biológicas?
Para poder hacer el análisis de expresión diferencial, es necesario normalizar nuestros datos, para esto edgeR nos permite calcular un factor de normalización.
```r
dgeNorm = calcNormFactors(dge)

#Puedes visualizar los factores de normalización
dgeNorm$samples
```

<img width="935" height="285" alt="Captura de pantalla 2025-08-27 a la(s) 19 08 15" src="https://github.com/user-attachments/assets/fe28e911-1514-4736-ad08-a55211c6e9ef" />

## Realizando el análisis de expresión diferencial

```r
#estimaremos la dispersión entre nuestras librerías 
dgeNorm = estimateCommonDisp(dgeNorm)

#podemos visualizar este valor (si no tuvieran replicas es posible ingresar este valor manualmente, pero no es ideal)
dgeNorm$common.dispersion
```
<img width="929" height="50" alt="Captura de pantalla 2025-08-27 a la(s) 19 09 01" src="https://github.com/user-attachments/assets/1a920bbc-724a-400d-a282-5bceb5926596" />

## Calculando la expresión diferencial entre nuestras librerias (prueba pareada, exact-test)
```r
#Pueden especificar que grupo de datos desean comparar
diff_exp = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("wt_sc", "wt_sl" ))
diff_exp2 = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("st_sc", "st_sl" ))

#ahora visualizaremos el objeto resultante de la prueba exacta de Fisher (objeto DGEExact)

diff_exp
```
<img width="928" height="311" alt="Captura de pantalla 2025-08-27 a la(s) 19 09 53" src="https://github.com/user-attachments/assets/662c8e73-55e6-4fb0-8ad5-59833d49108c" />


```r
#si desean saber más sobre esta función
?exactTest

#podemos preguntar las dimenciones del objeto resultante
dim(diff_exp)

```

<img width="937" height="49" alt="Captura de pantalla 2025-08-27 a la(s) 19 10 29" src="https://github.com/user-attachments/assets/b8b769da-eb88-4176-a6a0-07bd07734282" />

```r
#es deseable saber que genes tienen un mayor cambio y con que valor

topTags(diff_exp)
```


<img width="933" height="257" alt="Captura de pantalla 2025-08-27 a la(s) 19 12 31" src="https://github.com/user-attachments/assets/ca93c736-a90b-4e03-bf10-43c816f328c3" />


```r
#podemos conocer cuantos genes nos esta arrojando el topTags

dim(topTags(diff_exp))

```

<img width="925" height="50" alt="Captura de pantalla 2025-08-27 a la(s) 19 19 04" src="https://github.com/user-attachments/assets/d106e809-4176-4944-97f7-ae3881dd5183" />


## Guardar la tabla con todos los valores en un objeto de R y poder usarlo posteriormente es importante para nuestros fines.

```r

deTab = topTags(diff_exp, n=Inf)$table

# esta tabla puede ser tratada con los comandos de selección y comparación que ya hemos visto

deTab[c(15,30),]

```


<img width="933" height="86" alt="Captura de pantalla 2025-08-27 a la(s) 19 15 06" src="https://github.com/user-attachments/assets/a4be8ca8-904d-408c-8f02-423fd50c0e1f" />

```r
#podemos filtrar esta tabla por cualquier criterio matemático

row.names(deTab)[deTab$logFC > 5]  #recuerden que el FC esta dado en log2.
```

<img width="937" height="123" alt="Captura de pantalla 2025-08-27 a la(s) 19 21 49" src="https://github.com/user-attachments/assets/960faedd-09a4-40ed-8386-37d75f2a8c95" />

```r
#y podemos explorar los valores 
deTab["YNL284C-A",]

```

