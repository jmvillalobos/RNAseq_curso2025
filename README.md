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

```r
tabla <- read.table("https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/main/Saccharomyces.txt",
                    header = TRUE, sep = "\t")

```

Para comenzar a trabajar con nuestros datos de expresión debemos posicionarnos en el directorio donde esta nuestra tabla de datos, para posteriormente cargar la tabla en R.

Llamaremos a nuestra tabla y comenzaremos a manipularla
