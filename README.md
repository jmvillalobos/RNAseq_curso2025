# RNAseq curso 2025
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
<img width="941" height="72" alt="Captura de pantalla 2025-08-27 a la(s) 19 23 56" src="https://github.com/user-attachments/assets/b3d79611-0b74-4360-94ee-64ebef162a96" />



```r
#esto nos permite entonces filtrar automáticamente por cualquier valor que predeterminemos
#aquí se proponen valores estándar de filtrado pero esto se puede modificar

deGenes = rownames(deTab)[deTab$FDR < 0.05 & abs(deTab$logFC) > 1]
down=row.names(deTab)[deTab$logFC< -1]
over=row.names(deTab)[deTab$logFC> 1]

#para saber el número total de genes que pasan estos filtros (genes diferenciales según nuestro criterio)

print(paste("total de diferenciales:", length(deGenes)))

```



<img width="938" height="248" alt="Captura de pantalla 2025-08-27 a la(s) 19 24 38" src="https://github.com/user-attachments/assets/a3176107-858e-413f-ad41-fd1f93ddb311" />

## Graficando los genes diferenciales para poder observar cómo se comportan y guardando los resultados

```r
plotSmear(dge, de.tags=deGenes, ylab = "WT-sc vs WT-sl")
```

<img width="686" height="423" alt="Captura de pantalla 2025-08-27 a la(s) 19 26 32" src="https://github.com/user-attachments/assets/c8238b89-d4e8-40bc-9618-db1348c11ccb" />

```r

#Podemos exportar la tabla de genes diferenciales a nuestra carpeta de trabajo

write.table(deTab, file=paste(outpath, "diff_gene_wt.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

#aquí exportamos solo los genes que se sobreexpresan o se reprimen.
write.table(down, file=paste(outpath, "down.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(over, file=paste(outpath, "up.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

```


# Gráficos de expresión diferencial
En esta sección les mostraré algunos de los gráficos con los que regularmente se presentan los resultados de los análisis que hemos realizado. Son básicos y ustedes pueden encontrar en la red ejemplos más sofisticados, todo depende de que quieran mostrar. Las paqueterías de R que usaremos serán gplots.

```r
#Instalando las paqueterias
install.packages("gplots")
install.packages("RcolorBrewer")

library("gplots")
library("RColorBrewer")

```

# Visualizando nuestros datos de expresión diferencial en otro formato más interesante

```r
#normalizamos nuestros datos de expresión por cuentas por millón
normalizados= cpm(counts)

#extraemos la expresión de los genes diferenciales
normalizados_diferenciales= normalizados[deGenes,]

#veamos cómo se ve esta tablita
head(normalizados_diferenciales)

```
<img width="943" height="295" alt="Captura de pantalla 2025-08-27 a la(s) 19 31 40" src="https://github.com/user-attachments/assets/5f993bc7-892c-44a8-94c9-626d11f23fb2" />


```r
#hagamos un HEATMAP!
heatmap(normalizados_diferenciales)

```

<img width="495" height="481" alt="Captura de pantalla 2025-08-27 a la(s) 19 32 40" src="https://github.com/user-attachments/assets/a5bc9253-17cd-40fe-997d-f2654b493756" />


# Comparación de normalizados contra no-normalizados

```r
par(mfrow=c(1,2)) 

boxplot(log(counts),col=rainbow(6), main="antes de la normalización")
boxplot(log(normalizados), col=rainbow(6), main="después de la normalización")

```

<img width="659" height="452" alt="Captura de pantalla 2025-08-27 a la(s) 19 33 56" src="https://github.com/user-attachments/assets/1b55170e-21fd-4b4c-859c-6488d8265c9b" />

*En estos gráficos vemos algunas alarmas, dado que tenemos algunos ceros que al ser convertidos por log2, se van al infinito.*


# Barplot de la cantidad de reads en cada librería

```r

barplot(apply(normalizados_diferenciales,2,sum),las=2, cex.names = 1, col = (1:6))

```

<img width="697" height="430" alt="Captura de pantalla 2025-08-27 a la(s) 19 36 40" src="https://github.com/user-attachments/assets/923e75bf-6628-4d3a-972d-8755b5955c07" />

# Un PCA bonito de la distribución de nuestros datos

```r

pca <- princomp(normalizados_diferenciales[,c(1:6)])
plot(pca$loadings, col=as.factor(colnames(normalizados_diferenciales[,c(1:6)])),  pch=19, cex=2, main="con nitrógeno")
text(pca$loadings, as.vector(colnames(normalizados_diferenciales[,c(1:6)])), pos=3, cex=0.8)

```

<img width="665" height="470" alt="Captura de pantalla 2025-08-27 a la(s) 19 37 35" src="https://github.com/user-attachments/assets/acd3f5f2-afa6-4ae6-8652-a1572630dd8b" />


```r
pca <- princomp(normalizados_diferenciales[,c(7:12)])
plot(pca$loadings, col=as.factor(colnames(normalizados_diferenciales[,c(7:12)])),  pch=19, cex=2, main="si nitrógeno")
text(pca$loadings, as.vector(colnames(normalizados_diferenciales[,c(7:12)])), pos=3, cex=0.8)

```

<img width="709" height="452" alt="Captura de pantalla 2025-08-27 a la(s) 19 38 07" src="https://github.com/user-attachments/assets/0f54f51a-4f11-41b0-b735-6852fd8a908a" />

# Hagamos unos bonitos volcano-plot, nos ayudaran a comparar los datos

```r

with(deTab, plot(logFC, -log10(FDR), pch=20, cex=0.8, col="black", main="WT+N vs WT-N", xlim=c(-8, 8), ylim=c(0,300)))
text(deTab[1:20,]$logFC,-log(deTab[1:20,]$FDR,10),labels=rownames(deTab[1:20,]),cex=0.7,pos=1)
with(subset(deTab, FDR<.01 & abs(logFC)>2), points(logFC, -log10(FDR), pch=20, cex=0.5, col="green"))
abline(v=2,lty=2, col="blue")
abline(v=-2,lty=2, col="blue")
legend("bottomright","Up_regulated",cex=1)
legend("bottomleft","Down_regulated",cex=1)

```


<img width="756" height="456" alt="Captura de pantalla 2025-08-27 a la(s) 19 39 27" src="https://github.com/user-attachments/assets/a9784b64-b328-4313-abaa-020feddcf6e5" />

```r
deTab2 = topTags(diff_exp2, n=Inf)$table
deGenes2 = rownames(deTab2)[deTab2$FDR < 0.05 & abs(deTab2$logFC) > 1]


with(deTab2, plot(logFC, -log10(FDR), pch=20, cex=0.8, col="black", main="ste12+N vs ste12-N", xlim=c(-10, 10), ylim=c(0,320)))
text(deTab2[1:20,]$logFC,-log(deTab2[1:20,]$FDR,10),labels=rownames(deTab2[1:20,]),cex=0.7,pos=1)
with(subset(deTab2, FDR<.01 & abs(logFC)>2), points(logFC, -log10(FDR), pch=20, cex=0.5, col="green"))
abline(v=2,lty=2, col="blue")
abline(v=-2,lty=2, col="blue")
legend("bottomright","Up_regulated",cex=1)
legend("bottomleft","Down_regulated",cex=1)
```

<img width="722" height="469" alt="Captura de pantalla 2025-08-27 a la(s) 19 40 24" src="https://github.com/user-attachments/assets/6791983c-2947-45ac-af5a-d01db8ca2ec5" />


Si desean saber qué genes no se comparten entre ambas cepas en los inducidos, podemos hacer lo siguiente

```r
WTover= head(rownames(deTab), 30)
ste12over= head(rownames(deTab2), 30)

setdiff(WTover, ste12over)

```

<img width="930" height="53" alt="Captura de pantalla 2025-08-27 a la(s) 19 41 10" src="https://github.com/user-attachments/assets/eaddcb03-4d39-4d7a-a997-dbad70b3a3a3" />


## Ahora veremos algunas herramientas en línea que nos pueden ayudar a darle más sentido a nuestro análisis:

https://genemania.org/

https://www.yeastgenome.org/cgi-bin/GO/goTermFinder.pl

# Súper!!!! ya saben hacer análisis de expresión diferencial

