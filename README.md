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


## An object of class "DGEList"
## $counts
##         wt_sc_1 wt_sc_2 wt_sc_3 st_sc_1  st_sc_2 st_sc_3 wt_sl_1 wt_sl_2
## YPR161C    1061    1375    1321    1189 1285.000    1334    1158    1409
## YOL138C     552     694     694     589  616.072     671    1581    1871
## YDR395W    2691    3331    3741    3310 3281.000    3473    1403    1649
## YGR129W     352     513     350     444  418.000     319     339     388
## YPR165W    2862    3900    3101    3749 3763.000    3274    3377    3899
##         wt_sl_3 st_sl_1 st_sl_2 st_sl_3
## YPR161C    1070  1300.0  1379.0    1444
## YOL138C    1411  1665.4  1619.2    1951
## YDR395W    1297  1471.0  1651.0    1975
## YGR129W     266   381.0   334.0     424
## YPR165W    2926  3797.0  3533.0    3825
## 5826 more rows ...
## 
## $samples
##         group lib.size norm.factors
## wt_sc_1 wt_sc  9934217            1
## wt_sc_2 wt_sc 12547058            1
## wt_sc_3 wt_sc 11625013            1
## st_sc_1 st_sc 12644400            1
## st_sc_2 st_sc 12841546            1
## 7 more rows ...




