Se crea una lista para identificar los quants obtenidos en la cuantificacion realizada mediante salmon
```
samples <- list.files(path = "./quants", full.names = T, pattern="ERR") 
samples
```
Instalacion de los paquetes necesarios
```
install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("stringr")
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2,DEGreport,tximport")
```
Librerias a utilizar en el procesamiento 
```
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tibble)
library(apeglm)
```
Extrae la informacion de las cuantificaciones de cada una de las muestras 
```files <- file.path(samples, "quant.sf")
```
files
Se reemplaza en los nombres de las muestras quitando la terminacion quants
```
names(files) <- str_replace(samples, "./quants/", "") 
names(files) <- str_replace(names(files), "_quant", "") 
files
```
Se lee el archivo indice de los genes de homo sapiens 
```
tx2gene <- read.delim("tx2gene_grch38_ens94.txt")
```
Se observa la informacion de los genes presentes  
```
tx2gene %>% View()
```
Se aplica tximport para concatenar nuestros datos de las muestras con los genes de homo sapiens
```
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion=TRUE)
file.exists(files)
attributes(txi) #Se muestran los atributos del conjunto creado mediante tximport
```
```
Verificar si hay valores NA
anyNA(txi$counts)
Verificar si hay valores NaN
any(is.nan(txi$counts))
```
Se guarda los datos y los genes en un objeto tipo data frame
```
data <- txi$counts %>% 
  round(digits = 2) %>%  # Redondear a dos decimales
  as.data.frame()
```
Se importa los datos de la metadata
```
samples <- read.table("C:/Users/Det-Pc/Documents/Maestria UNIR/TFM/RStudio - Results/SraRunTable.csv", header=TRUE, sep=",", quote="")
```
Se verifica que las filas y columnas son del mismo tamaño
```
ncol(txi$counts)
nrow(samples)
```
Se convierte a la columna growth del data frame a objeto tipo factor
```
samples$growth <- as.factor(samples$growth)
```
Se crea el objeto a clase deseqdatasetfromtximport, para llevar a cabo el analisis de expresion diferencial 
```
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ growth)
ddsTxi
```
Se toma como linea base las muestras de normal tissue 
```
ddsTxi$growth <- relevel(ddsTxi$growth, ref = "normal tissue")
```
Resumen de las estadisticas descriptivas de la base de datos cuantificados
``
summary(counts(ddsTxi))
```
Se aplica el analisis de expresion diferencial
```
ddsTxi <- DESeq(ddsTxi)
res <- results(ddsTxi)
res
resdata <- as.data.frame(res)
resdatagene <- rownames(resdata)
```
Se guarda en un archivo los resultados obtenidos
```
write.csv(resdatagene, file = "genes.csv", row.names = FALSE)
```
Debido a la reutilizacion de los resultados se los guarda en otras variables. 
```
resvolcano <- res
ddstxi2 <- ddsTxi #Se crea en la copia del data frame de los datos 
```
Agrupar los datos por sexo y contar los individuos únicos para generar las graficas que describen la poblacion de estudio
```
count_individuals <- samples %>%
  group_by(sex) %>%
  summarise(count = n_distinct(Individual))
Definir una paleta de colores
colores1 <- c("#3d8b7d", "#8fbc91", "#dbc557")
Crear el gráfico con ggplot2
ggplot(count_individuals, aes(x = sex, y = count, fill = sex)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5) +
  scale_fill_manual(values = colores1) +
  labs(title = "Número de Individuos por Sexo", x = "Sexo", y = "Número de Individuos") +
  theme_minimal()
#Agrupar los datos por AGE y contar los individuos únicos
count_individuals <- samples %>%
  group_by(AGE) %>%
  summarise(count = n_distinct(Individual)) %>%
  mutate(AGE = as.factor(AGE))  # Convertir AGE a variable categórica
#Definir una paleta de colores
colores1 <- c("#2a9d8f", "#e9c46a", "#f4a261", "#264653", "#f1faee", "#3aafa9", "#72efdd", 
              "#e07a5f", "#2a4d69", "#f1faee", "#c5d0e3")
#Crear el gráfico con ggplot2
ggplot(count_individuals, aes(x = AGE, y = count, fill = AGE)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5) +
  scale_fill_manual(values = colores1) +
  labs(title = "Número de Individuos por Edad", x = "Edad", y = "Número de Individuos") +
  theme_minimal()
#Calcular el conteo y el porcentaje de cada tipo de crecimiento por tipo de muestra
#Crear el gráfico de barras apiladas horizontal con el número de muestras
ggplot(samples, aes(x = sampletype, fill = growth)) +
  geom_bar(position = "stack") +
  geom_text(stat = 'count', aes(label = ..count..), position = position_stack(vjust = 0.5), color = "white") + scale_fill_manual(values = colores1) +
  labs(title = "Distribución de Condición de Crecimiento y el Sitio de muestreo",
       y = "Número de Muestras",
       x = "Sitio de Muestreo") +
  theme_minimal()
```
Normalizacion de datos
ddstxi2 <- estimateSizeFactors(ddstxi2)
#Normalizacion mediante logaritmo regularizado
dds_rlog <- rlog(ddstxi2)
#Normalizacion mediante transformacion estabilizadora de varianza 
dds_vst <- varianceStabilizingTransformation(ddstxi2)
assay(ddstxi2, "rlog") <- assay(dds_rlog)
assay (ddstxi2, "vst") <- assay(dds_vst)
#Se guarda los datos y los datos normalizados en formato .rds
saveRDS(ddstxi2, "ddstxi2.rds", compress = T) 
saveRDS(dds_rlog, "dds_rlog.rds", compress = T)
saveRDS(dds_vst, "dds_vst.rds", compress=T)
saveRDS(samples, "samples.rds", compress=T)
#Creacion de carpeta para guardar los graficos a generar
if (!file.exists("QCplots")){
dir.create("QCplots")
}





#PCA pLots
# make my customized function to plot PCA
plotPCA_vst <- function (object, assay, ntop = 1000) {
  rv <- rowVars(assay(object, assay))
  select <- order(rv, decreasing = TRUE) [seq_len(min(ntop, length(rv) ))]
  pca <- prcomp(t(assay(object, assay) [select, ]), center=TRUE , scale=FALSE)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  df <- cbind( as.data.frame(colData(object)), pca$x)
  # order points so extreme samples are more Likely to get lab
  ord <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
  df <- df[ord,]
  attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar),percentVar=100*percentVar)
  return (df)
}
pca.data <- plotPCA_vst(ddstxi2, assay="vst", ntop=1000)
saveRDS(pca.data, "pca.data.rds", compress=T)
for (i in c(1:4)) {
  for (j in c( (i+1):5)) {
    percentVar <- round(attr(pca.data, "percentVar")$percentVar)
    p <- ggplot(pca.data, aes_string(x=paste0("PC",i), y=paste0("PC",j), color="growth")) +
      geom_point (size=3) + xlab(paste0("PC",i,": ",percentVar[i],"% variance")) +
      ylab(paste0("PC",j,": ",percentVar[j],"% variance")) +
      
      theme_bw()
    ggsave(filename = paste0("QCplots/PCAplot.","PC",i,"vs","PC",j,".pdf"), plot = p, width=6, height=5)
  }}





#Boxplot
#Prepara los datos normalizados 
counts_vst <- ddstxi2@assays@data$vst
dat <- counts_vst %>%
as.data.frame() %>%
rownames_to_column (var="gene") %>%
gather (key = "Run", value = "counts_vst", -gene) %>%
left_join(samples, by="Run")
#Se crea una combinacion en los factores entre growth y run
dat <- dat %>%
  mutate(Run = factor(Run, levels = unique(Run[order(growth)])))
#Graficos de boxplot analizando las condiciones de la muestras en normalizacion VST 
p_1 <- ggplot(data=dat, aes(x=Run,y=counts_vst, fill=growth)) +
geom_boxplot() + labs (x="Muestras",y="Conteo de lecturas normalizadas VST", fill="Condición")+
ggtitle("Conteo de lecturas normalizadas VST por Muestra") +
theme (axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1))
p_1
ggsave(filename = "QCplots/boxplot-counts_vst.jpg", plot = p_1, width = 10, height = 5) #Se guarda el grafico en la carpeta creada
#Se prepara los datos para el gtrafico usando los datos normalizados mediante rlog
counts_rlog <- ddstxi2@assays@data$rlog
dat2 <- counts_rlog %>%
as.data.frame() %>%
rownames_to_column (var="gene") %>%
gather (key = "Run", value = "counts_rlog", -gene) %>%
left_join(samples, by="Run")
dat2 <- dat2 %>%
  mutate(Run = factor(Run, levels = unique(Run[order(growth)])))
#Genera el grafico de boxplot mediante los datos normalizados con rlog
p_2 <- ggplot(data=dat2, aes (x=Run,y=counts_rlog, fill=growth)) +
geom_boxplot() +
labs(x="Muestras" ,y="Conteo de lecturas normalizadas rlog", fill="Condición")+
ggtitle("Conteo de lecturas normalizadas rlog por Muestra") +
theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
p_2
ggsave(filename = "QCplots/boxplot-counts_rlog.jpg", plot = p_2, width = 10, height = 5)
#Generacion de los graficos tipo heatmap para comparacion entre muestras 
sampleDists <- dist(t(counts_vst)) #se usa los datos normalziados mediante vst
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Oranges")) )(255) #Paleta de colores para el grafico 
mydata_col = data.frame(Condición=samples$growth, row.names = samples$Run)
latevsnormal <- mydata_col %>% filter(Condición %in% c("normal tissue", "late organoid culture"))
earlyvsnormal <- mydata_col %>% filter(Condición %in% c("normal tissue", "early organoid culture"))
neoplasmavsnormal <- mydata_col %>% filter(Condición %in% c("normal tissue", "neoplasm"))
latevsnormal <- data.frame(latevsnormal)
#Se genera el grafico pheatmap
p_3 <- pheatmap(
sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors, 
annotation_col = mydata_col,
main="Distancia Euclidiana entre muestras"
)
p_3
#Se guarda el grafico en la carpeta creada anteriormente 
ggsave(filename = "QCplots/sample_heatmap.jpg", plot = p_3, width = 10, height = 8)
#Se crea una copia de los datos 
ddstxi3 <- ddstxi2
resultsNames(ddstxi3)
#Se toma los valores de normalizacion para generar los histogramas   
LFC1 <- lfcShrink(ddstxi3, coef = "growth_late.organoid.culture_vs_normal.tissue", type = "apeglm")
LFC2 <- lfcShrink(ddstxi3, coef = "growth_early.organoid.culture_vs_normal.tissue", type = "apeglm")
LFC3 <- lfcShrink(ddstxi3, coef = "growth_neoplasm_vs_normal.tissue", type = "apeglm")

subexpresados1 <- LFC1[!is.na(LFC1$padj) & LFC1$padj < 0.05 & LFC1$log2FoldChange < -2, ]
genes_subexpresados1 <- data.frame(GeneName = rownames(subexpresados1),
                                  log2FoldChange = subexpresados1$log2FoldChange,
                                  padj = subexpresados1$padj)
subexpresados2 <- LFC2[!is.na(LFC2$padj) & LFC2$padj < 0.05 & LFC2$log2FoldChange < -2, ]
genes_subexpresados2 <- data.frame(GeneName = rownames(subexpresados2),
                                  log2FoldChange = subexpresados2$log2FoldChange,
                                  padj = subexpresados2$padj)
subexpresados3 <- LFC3[!is.na(LFC3$padj) & LFC3$padj < 0.05 & LFC3$log2FoldChange < -2, ]
genes_subexpresados3 <- data.frame(GeneName = rownames(subexpresados3),
                                  log2FoldChange = subexpresados3$log2FoldChange,
                                  padj = subexpresados3$padj)

sobreexpresados1 <- LFC1[!is.na(LFC1$padj) & LFC1$padj < 0.05 & LFC1$log2FoldChange > 2, ]
genes_sobreexpresados1 <- data.frame(GeneName = rownames(sobreexpresados1),
                                   log2FoldChange = sobreexpresados1$log2FoldChange,
                                   padj = sobreexpresados1$padj)
sobreexpresados2 <- LFC2[!is.na(LFC2$padj) & LFC2$padj < 0.05 & LFC2$log2FoldChange > 2, ]
genes_sobreexpresados2 <- data.frame(GeneName = rownames(sobreexpresados2),
                                     log2FoldChange = sobreexpresados2$log2FoldChange,
                                     padj = sobreexpresados2$padj)
sobreexpresados3 <- LFC3[!is.na(LFC3$padj) & LFC3$padj < 0.05 & LFC3$log2FoldChange > 2, ]
genes_sobreexpresados3 <- data.frame(GeneName = rownames(sobreexpresados3),
                                     log2FoldChange = sobreexpresados3$log2FoldChange,
                                     padj = sobreexpresados3$padj)
genes1 <- genes_sobreexpresados1
genes2 <- genes_sobreexpresados2
genes3 <- genes_sobreexpresados3
genes4 <- genes_subexpresados1
genes5 <- genes_subexpresados2
genes6 <- genes_subexpresados3
genes1$source <- "lfc1"
genes2$source <- "lfc2"
genes3$source <- "lfc3"
genes4$source <- "lfc1"
genes5$source <- "lfc2"
genes6$source <- "lfc3"
# Combinar los data frames
genessobreexpresados <- rbind(genes1, genes2, genes3)
genessubreexpresados <- rbind(genes4, genes5, genes6)

write.csv(genessobreexpresados, file = "genessobreexpresados.csv", row.names = FALSE)
write.csv(genessubreexpresados, file = "genessubreexpresados.csv", row.names = FALSE)
#Crear y guardar el histograma usando funciones de gráficos base de R
jpeg(filename = "QCplots/latevsnormal.jpg", width = 12, height = 8, units = "in", res = 300)
hist(LFC1$pvalue, breaks = 50, col = 'grey', 
     main = NULL, 
     xlab = 'p-value',
     cex.lab = 1.5, 
     cex.axis = 1.5)
dev.off()
#Histograma realizado con los valores p ajustados 
jpeg(filename = "QCplots/latevsnormaladj.jpg", width = 10, height = 6, units = "in", res = 300)
hist(LFC1$padj, breaks = 50, col = 'grey', 
     main = NULL, 
     xlab = 'Adjusted p-value',
     cex.lab = 1.5, 
     cex.axis = 1.5)
dev.off()
#Se puede verificar cual es el numero de genes sobreexpresados o subexpresados
attach(as.data.frame(LFC1))
#Numero total de genes ajustados con p-value<0.05
summary(LFC1, alpha=0.05)
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
#Subexpresion
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) 
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) 
#Sobreexpresion
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) 
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) 
#Se crea el grafico histograma para organoide temprano vs tejido normal 
jpeg(filename = "QCplots/tempranovsnormal.jpg", width = 10, height = 6, units = "in", res = 300)
hist(LFC2$pvalue, breaks = 50, col = 'grey', 
     main = NULL, 
     xlab = 'p-value',
     cex.lab = 1.5, 
     cex.axis = 1.5)
dev.off()
#Histograma con los valores p ajustados
jpeg(filename = "QCplots/tempranovsnormaladj.jpg", width = 10, height = 6, units = "in", res = 300)
hist(LFC2$padj, breaks = 50, col = 'grey', 
     main = NULL, 
     xlab = 'Adjusted p-value',
     cex.lab = 1.5, 
     cex.axis = 1.5)
dev.off()
#Calculo de genes expresados para LFC2
attach(as.data.frame(LFC2))
#Numero de genes expresados con valores ajustados p-value<0.05
summary(LFC2, alpha=0.05)
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
#Subexpresados
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) 
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) 
#Sobreexpresion
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) 
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) 

#Distribucion de p-values de LCF3
jpeg(filename = "QCplots/neoplasmavsnormal.jpg", width = 10, height = 6, units = "in", res = 300)
hist(LFC3$pvalue, breaks = 50, col = 'grey', 
     main = NULL, 
     xlab = 'p-value',
     cex.lab = 1.5, 
     cex.axis = 1.5)
dev.off()
#Histograma mediante datos de valores p ajustados
jpeg(filename = "QCplots/neoplasmavsnormaladj.jpg", width = 10, height = 6, units = "in", res = 300)
hist(LFC3$padj, breaks = 50, col = 'grey', 
     main = NULL, 
     xlab = 'Adjusted p-value',
     cex.lab = 1.5, 
     cex.axis = 1.5)
dev.off()
#Genes expresados diferencialmente
attach(as.data.frame(LFC3))
summary(LFC3, alpha=0.05)
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
#Subexpresion
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) 
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) 
#Sobrexpresion
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) 
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) 
#Convertir los resultados a data frames y agregar una columna que indique el contraste
LFC1_df <- as.data.frame(LFC1)
LFC2_df <- as.data.frame(LFC2)
LFC3_df <- as.data.frame(LFC3)
# Identificar los genes con una expresión diferencial significativa
DEGs1 <- LFC1_df[!is.na(LFC1_df$pvalue) & LFC1_df$pvalue < 0.05 & abs(LFC1_df$log2FoldChange) > 2, ]
DEGs2 <- LFC2_df[!is.na(LFC2_df$pvalue) & LFC2_df$pvalue < 0.05 & abs(LFC2_df$log2FoldChange) > 2, ]
DEGs3 <- LFC3_df[!is.na(LFC3_df$pvalue) & LFC3_df$pvalue < 0.05 & abs(LFC3_df$log2FoldChange) > 2, ]
# Separar genes con aumento y disminución de expresión
increased_genes1 <- rownames(DEGs1[DEGs1$log2FoldChange > 2, ])
decreased_genes1 <- rownames(DEGs1[DEGs1$log2FoldChange < -2, ])
increased_genes2 <- rownames(DEGs2[DEGs2$log2FoldChange > 2, ])
decreased_genes2 <- rownames(DEGs2[DEGs2$log2FoldChange < -2, ])
increased_genes3 <- rownames(DEGs3[DEGs3$log2FoldChange > 2, ])
decreased_genes3 <- rownames(DEGs3[DEGs3$log2FoldChange < -2, ])
# Combinar las listas de genes
representative_genes1 <- unique(c(increased_genes1, decreased_genes1))
representative_genes2 <- unique(c(increased_genes2, decreased_genes2))
representative_genes3 <- unique(c(increased_genes3, decreased_genes3))
#Subconjunto de counts_vst para los genes representativos
desired_columns1 <- rownames(latevsnormal)  # Esto toma los nombres de las columnas 
counts_vst_representative1 <- counts_vst[rownames(counts_vst) %in% representative_genes1, colnames(counts_vst) %in% desired_columns1]
desired_columns2 <- rownames(earlyvsnormal)  # Esto toma los nombres de las columnas 
counts_vst_representative2 <- counts_vst[rownames(counts_vst) %in% representative_genes2, colnames(counts_vst) %in% desired_columns2]
desired_columns3 <- rownames(neoplasmavsnormal)  # Esto toma los nombres de las columnas 
counts_vst_representative3 <- counts_vst[rownames(counts_vst) %in% representative_genes3, colnames(counts_vst) %in% desired_columns3]


as.matrix()
slice(1:1000) %>%
  counts_vst_var <- counts_vst %>%
  counts_vst_var <- counts_vst[!is.na(counts_vst$) & counts_vst$pvalue < 0.05 & abs(counts_vst$log2FoldChange) > 2, ]
p_4 <- pheatmap(counts_vst_var [,1:35], scale="row", color = colorRampPalette(c("darkblue", "white","darkred")) (1000), annotation_col = mydata_col,main="Gene Expression (VST) of Top 1000 Variable Genes",show_rownames = F)


#Crear un mapa de calor
p_6 <- pheatmap(counts_vst_representative1,
         scale = "row",
         color = colorRampPalette(c("darkblue", "white", "darkred"))(50),
         annotation_col = latevsnormal,
         main = "Gene Expression of Representative DEGs",
         show_rownames = TRUE)
ggsave(filename = "QCplots/latevsnormalheat.jpg", plot = p_6, width = 10, height = 10)
p_7 <- pheatmap(counts_vst_representative2,
                scale = "row",
                color = colorRampPalette(c("darkblue", "white", "darkred"))(50),
                annotation_col = earlyvsnormal,
                main = "Gene Expression of Representative DEGs",
                show_rownames = TRUE)
ggsave(filename = "QCplots/earlyvsnormalheat.jpg", plot = p_7, width = 10, height = 12)
p_8 <- pheatmap(counts_vst_representative3,
                scale = "row",
                color = colorRampPalette(c("darkblue", "white", "darkred"))(50),
                annotation_col = neoplasmavsnormal,
                main = "Gene Expression of Representative DEGs",
                show_rownames = TRUE)
ggsave(filename = "QCplots/neoplasmavsnormalheat.jpg", plot = p_8, width = 10, height = 10)

#KEGG pathways
BiocManager::install(c("pathview", "gage", "gageData"))
library("AnnotationDbi")
library("org.Hs.eg.db")
library(pathview)
library(gage)
library(gageData)
#Descargar la informacion de la libreria
columns(org.Hs.eg.db)

#Añadir a la matriz los datos de symbol, entrez y names de los genes
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

head(res, 10)

##descargar data
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
#determinar valores propios de los datos analizados
foldchanges = res$log2FoldChange
names(foldchanges) =res$entrez
head(foldchanges)
# Obtener los resultados de Kegg
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
#resumen de las vias reguladas positiva y negativamente
lapply(keggres, head)
# Obtener las vias principales
# Seleccionar los 5 pathways más significativos sobreexpresados
keggrespathways_greater = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  as_tibble() %>% 
  filter(row_number() <= 5) %>% 
  .$id %>% 
  as.character()
# Seleccionar los 5 pathways más significativos subexpresados
keggrespathways_less = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  as_tibble() %>% 
  filter(row_number() <= 5) %>% 
  .$id %>% 
  as.character()

# Combinar ambos conjuntos de pathways
keggrespathways = c(keggrespathways_greater, keggrespathways_less)

# Ver los pathways seleccionados
keggrespathways

# IDs de las vias
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
# Definir la funcion para graficar las rutas
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Funcion para graficar multiples graficas que esten en el vector keggresids
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

##
####ENRICHMENT
BiocManager::install("clusterProfiler")
library(clusterProfiler)
# Filtrar los genes significativos (ya hecho)
dsttgenes <- subset(res, (log2FoldChange > 2 | log2FoldChange < -2) & padj <= 0.05)
# Verifica qué tipo de ID estás usando en dsttgenes
head(rownames(dsttgenes))
write.csv(dsttgenes, file = "goenrichgenes.csv", row.names = FALSE)
# Realiza el análisis de enriquecimiento GO usando el keyType adecuado
GO_enrich <- enrichGO(gene = rownames(dsttgenes), # IDs de los genes (por ejemplo, ENSEMBL IDs)
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENSEMBL",        # Cambia a SYMBOL, ENTREZID, etc. según corresponda
                      ont = "MF",                 # Ontología de interés
                      minGSSize = 30)
write.csv(GO_enrich, file = "goenrich.csv", row.names = FALSE)
# Ver los resultados
summary(GO_enrich)

View(GO_enrich@result)

# Configurar los márgenes para que las etiquetas se vean bien
par(mar=c(5, 10, 3, 3))

# Crear el gráfico de barras
jpeg(filename = "QCplots/barplotgo.jpg", width = 10, height = 10, units = "cm", res = 300)
barplot(rev(-log10(GO_enrich@result$p.adjust[1:10])),     # -log10 de los p-values ajustados
        horiz = TRUE,                                     # Barras horizontales
        names=rev(GO_enrich@result$Description[1:10]),    # Descripción de los términos GO
        las=2,                                            # Rotar las etiquetas de los ejes
        xlab="-log10(adj.p-value)",                       # Etiqueta del eje x
        cex.names = 0.7,                                  # Tamaño de las etiquetas
        col="#dbc557")                                 # Color de las barras
dev.off()



# Ajustar márgenes y aumentar el tamaño de la imagen
jpeg(filename = "QCplots/barplotgo.jpg", width = 20, height = 8, units = "in", res = 300)

# Ajustar los márgenes para evitar que las etiquetas se corten
par(mar = c(5, 50, 4, 2) + 0.1, cex.lab = 1.5)  # Márgenes: abajo, izquierda, arriba, derecha

# Crear el barplot
barplot(rev(-log10(GO_enrich@result$p.adjust[1:10])),    # -log10 de los p-values ajustados
        horiz = TRUE,                                    # Barras horizontales
        names=rev(GO_enrich@result$Description[1:10]),   # Descripción de los términos GO
        las=2,                                           # Rotar las etiquetas de los ejes
        xlab="-log10(adj.p-value)",                      # Etiqueta del eje x
        cex.names = 2,
        cex.axis = 1.5,                                  # Tamaño de las etiquetas
        col="#dbc557")                                   # Color de las barras
# Añadir una línea vertical en el -log10(0.05) como referencia
abline(v=-log10(0.05))
dev.off()

# Dotplot on enrichResult and gseaResult objects:
jpeg(filename = "QCplots/goenrich.jpg", width = 6, height = 4, units = "in", res = 300)
enrichplot::dotplot(GO_enrich, orderBy="p.adjust")
dev.off()
BiocManager::install("EnhancedVolcano", force = TRUE)
library(EnhancedVolcano)
  
  LFC1$symbol <- tx2gene$symbol[match(row.names(LFC1), tx2gene$ensgene)]
  LFC2$symbol <- tx2gene$symbol[match(row.names(LFC2), tx2gene$ensgene)]
  LFC3$symbol <- tx2gene$symbol[match(row.names(LFC3), tx2gene$ensgene)]
  vol_lab1 <- EnhancedVolcano(LFC1,
                             lab = LFC1$symbol,
                             labSize = 4,
                             x= 'log2FoldChange',
                             y = 'padj',
                             ylab = bquote(~-Log [10] ~ italic(padj)),
                             pCutoff = 0.05,
                             xlim = c(-7,7),
                             ylim = c(0,5),
                             FCcutoff = 2,
                             title = paste0("Volcano plot"),
                             legendLabels=c('Genes no DE y cambio \nen la expresión menor a 2',
                                            'Genes no DE y cambio \nen la expresión mayor a 2',
                                            'Genes son DE y cambio \nen la expresión menor a 2',
                                            'Genes son DE y cambio \nen la expresión mayor a 2'),
                             legendPosition = 'right')
  vol_lab1
  ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-1",".png"), plot=vol_lab1, width = 14, height=7)      
  vol_lab2 <- EnhancedVolcano(LFC2,
                              lab = LFC2$symbol,
                              labSize = 4,
                              x= 'log2FoldChange',
                              y = 'padj',
                              ylab = bquote(~-Log [10] ~ italic(padj)),
                              pCutoff = 0.05,
                              xlim = c(-5,10),
                              ylim = c(0,6),
                              FCcutoff = 2,
                              title = paste0("Volcano plot"),
                              legendLabels=c('Genes no DE y cambio \nen la expresión menor a 2',
                                             'Genes no DE y cambio \nen la expresión mayor a 2',
                                             'Genes son DE y cambio \nen la expresión menor a 2',
                                             'Genes son DE y cambio \nen la expresión mayor a 2'),
                              legendPosition = 'right')
  vol_lab2
  ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-2",".png"), plot=vol_lab2, width = 14, height=7)  
  vol_lab3 <- EnhancedVolcano(LFC3,
                              lab = LFC3$symbol,
                              labSize = 4,
                              x= 'log2FoldChange',
                              y = 'padj',
                              ylab = bquote(~-Log [10] ~ italic(padj)),
                              pCutoff = 0.05,
                             # xlim = c(-10,10),
                              ylim = c(0,10),
                              FCcutoff = 2,
                              title = paste0("Volcano plot"),
                              legendLabels=c('Genes no DE y cambio \nen la expresión menor a 2',
                                             'Genes no DE y cambio \nen la expresión mayor a 2',
                                             'Genes son DE y cambio \nen la expresión menor a 2',
                                             'Genes son DE y cambio \nen la expresión mayor a 2'),
                              legendPosition = 'right')
vol_lab3
ggsave(filename = paste0("Volcanoplots/Volcanoplot_gene_name-3",".png"), plot=vol_lab3, width = 14, height=7)  
  
tmp = plotCounts(ddsTxi, gene = grep('ENSG00000163735', names(ddsTxi), value = TRUE), intgroup = "growth", pch = 18, main = '??? expression', returnData = TRUE)
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5))
p <- ggplot(tmp, aes(x = growth, y = count)) + geom_boxplot() + 
geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.6) + ggtitle('ENSG00000163735/ CXCL5 expression') +theme
print(p)

tmp1 = plotCounts(ddsTxi, gene = grep('ENSG00000230678', names(ddsTxi), value = TRUE), intgroup = "growth", pch = 18, main = '??? expression', returnData = TRUE)
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5))
p1 <- ggplot(tmp1, aes(x = growth, y = count)) + geom_boxplot() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.6) + ggtitle('ENSG00000230678/ BRD2 expression') +theme
print(p1)


tmp2 = plotCounts(ddsTxi, gene = grep('ENSG00000102038', names(ddsTxi), value = TRUE), intgroup = "growth", pch = 18, main = '??? expression', returnData = TRUE)
theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5))
p2 <- ggplot(tmp2, aes(x = growth, y = count)) + geom_boxplot() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.6) + ggtitle('ENSG00000102038/ SMARCA1 expression') +theme
print(p2)




# Instalar y cargar el paquete VennDiagram
install.packages("VennDiagram")
library(VennDiagram)

# Definir los conjuntos de datos
cultivo_tardio <- c("SLC7A2", "ABCB1", "CPVL", "AREG", "EREG", "FGF13", "GDF15", "THEMIS2", 
                    "AKAP12", "CDK7", "ALDH1L2", "SECTM1", "FHAD1", "CRABP2", "CREB5", "FSIP1", 
                    "PRDM8", "ADGRF1", "PLA2G4D", "CEL", "REG3A", "TNFRSF10C", "GJC1", "SLC44A4", 
                    "CSNK2B", "HLA-C", "MUC5AC", "EIF4BP7", "EIF3FP3", "HLA-B", "CFB", "MTND4P12", 
                    "PCDHGB6", "TRIL", "ALOX5", "AC256203.1", "CBWD4P", "CSNK2B", "ABHD16A")

cultivo_temprano <- c("SLC7A9", "EML1", "FBLN1", "ABCB1", "CERS4", "PALM", "GLRA2", "ZDHHC15", 
                      "CPVL", "B3GAT1", "SOX6", "LIN7A", "PROC", "XPNPEP2", "MAGEA10", "EREG", 
                      "SIX1", "HOXD10", "HOXD11", "GAMT", "SERPINF1", "ECHDC3", "ALDH1L2", "IGF2BP3", 
                      "NACAD", "CFAP300", "PLA2G12B", "CPLX2", "CREB5", "KIRREL3", "VENTX", "PDLIM3", 
                      "SAMD8", "SV2A", "AQP5", "ALOX15", "SLC1A7", "FREM1", "JSRP1", "FILIP1L", 
                      "APLN", "TCEA2", "ZNF556", "EGFL7", "CNTNAP2", "CD19", "PLAG1", "ZNF716", 
                      "DENND5A", "AC107956.1", "STMN3", "TPM2", "TGM2", "ONECUT3", "HLA-C", "LDHAP4", 
                      "RCN1P2", "COL28A1", "RPL13P12", "EIF4BP7", "MTND2P28", "AC098826.2", "BRD2", 
                      "EIF3FP3", "RPL10P9", "LDHAP7", "LDHAP3", "RPS4XP22", "MTND6P4", "ZKSCAN7", 
                      "MTRNR2L6", "FBXO17", "TAPBP", "OPRD1", "NOTCH1", "ADAMTS15")

neoplasma <- c("ZIC2", "ABCB11", "CACNG4", "FER1L4", "UPK3A", "WNT5B", "EREG", "CKMT2", "ECHDC3", 
               "TESPA1", "WDR17", "AKAP6", "AKR1C2", "CFAP221", "CXCL5", "CLDN2", "ARL4D", "EID2B", 
               "TSHZ2", "GJC1", "TMEM119", "PCDH9", "LSAMP", "HEATR4", "HLA-DQA1", "STMN3", "DUSP27", 
               "CES1", "RNA5SP199", "CR1", "CFB", "C2", "SLC44A4", "CSNK2B", "BAG6", "LDHAP4", 
               "EIF4BP7", "MTND2P28", "AC016734.1", "EIF3FP3", "LDHAP7", "LDHAP3", "RPP21", 
               "MTND4P12", "AC105052.3", "ZNF660", "RORA", "HLA-A")

# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    "Cultivo tardío" = cultivo_tardio,
    "Cultivo temprano" = cultivo_temprano,
    "Neoplasma" = neoplasma
  ),
  filename = "diagrama_venn.png",
  col = "transparent",
  fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.col = c("#E41A1C", "#377EB8", "#4DAF4A")
)

# Guardar el archivo de imagen
jpeg("diagrama_venn.jpg")
grid.draw(venn.plot)
dev.off()
# Instalar y cargar el paquete ggvenn
install.packages("ggvenn")
library(ggvenn)

# Definir los conjuntos de datos
cultivo_tardio <- c("SLC7A2", "ABCB1", "CPVL", "AREG", "EREG", "FGF13", "GDF15", "THEMIS2", 
                    "AKAP12", "CDK7", "ALDH1L2", "SECTM1", "FHAD1", "CRABP2", "CREB5", "FSIP1", 
                    "PRDM8", "ADGRF1", "PLA2G4D", "CEL", "REG3A", "TNFRSF10C", "GJC1", "SLC44A4", 
                    "CSNK2B", "HLA-C", "MUC5AC", "EIF4BP7", "EIF3FP3", "HLA-B", "CFB", "MTND4P12", 
                    "PCDHGB6", "TRIL", "ALOX5", "AC256203.1", "CBWD4P", "CSNK2B", "ABHD16A")

cultivo_temprano <- c("SLC7A9", "EML1", "FBLN1", "ABCB1", "CERS4", "PALM", "GLRA2", "ZDHHC15", 
                      "CPVL", "B3GAT1", "SOX6", "LIN7A", "PROC", "XPNPEP2", "MAGEA10", "EREG", 
                      "SIX1", "HOXD10", "HOXD11", "GAMT", "SERPINF1", "ECHDC3", "ALDH1L2", "IGF2BP3", 
                      "NACAD", "CFAP300", "PLA2G12B", "CPLX2", "CREB5", "KIRREL3", "VENTX", "PDLIM3", 
                      "SAMD8", "SV2A", "AQP5", "ALOX15", "SLC1A7", "FREM1", "JSRP1", "FILIP1L", 
                      "APLN", "TCEA2", "ZNF556", "EGFL7", "CNTNAP2", "CD19", "PLAG1", "ZNF716", 
                      "DENND5A", "AC107956.1", "STMN3", "TPM2", "TGM2", "ONECUT3", "HLA-C", "LDHAP4", 
                      "RCN1P2", "COL28A1", "RPL13P12", "EIF4BP7", "MTND2P28", "AC098826.2", "BRD2", 
                      "EIF3FP3", "RPL10P9", "LDHAP7", "LDHAP3", "RPS4XP22", "MTND6P4", "ZKSCAN7", 
                      "MTRNR2L6", "FBXO17", "TAPBP", "OPRD1", "NOTCH1", "ADAMTS15")

neoplasma <- c("ZIC2", "ABCB11", "CACNG4", "FER1L4", "UPK3A", "WNT5B", "EREG", "CKMT2", "ECHDC3", 
               "TESPA1", "WDR17", "AKAP6", "AKR1C2", "CFAP221", "CXCL5", "CLDN2", "ARL4D", "EID2B", 
               "TSHZ2", "GJC1", "TMEM119", "PCDH9", "LSAMP", "HEATR4", "HLA-DQA1", "STMN3", "DUSP27", 
               "CES1", "RNA5SP199", "CR1", "CFB", "C2", "SLC44A4", "CSNK2B", "BAG6", "LDHAP4", 
               "EIF4BP7", "MTND2P28", "AC016734.1", "EIF3FP3", "LDHAP7", "LDHAP3", "RPP21", 
               "MTND4P12", "AC105052.3", "ZNF660", "RORA", "HLA-A")

# Crear un diagrama de Venn con ggvenn
venn_data <- list("Cultivo tardío de organoide vs Tejido normal" = cultivo_tardio, 
                  "Cultivo temprano de organoide vs Tejido normal" = cultivo_temprano, 
                  "Neoplasma vs Tejido normal" = neoplasma)

# Crear el diagrama y ajustar las etiquetas
ggvenn(venn_data, show_elements = TRUE, fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"), text_size = 3, stroke_size = 1)

# Guardar la imagen
ggsave("diagrama_venn_con_genes.png", width = 10, height = 10, units = "in", dpi = 300)
