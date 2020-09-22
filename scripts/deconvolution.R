#!/usr/bin/env Rscript

# Need a single cell dataset -- gene expression rows vs cells and gene expression assigned to those cells. Then your bulk gene expression data

#-->1 Get top markers from single cell data
#-->2 of top markers which are best given some cutoff
#-->3 creating a reference matrix using those markers. What is the average expression for that marker gene for cells assigned to different cell types

############download Seurat

install.packages('Seurat')
library(Seurat)
library(Biobase)

##############get single cell and phenotype data 

cellmatrix = read.csv("single_cell_gene_matrix.csv")
rownames(cellmatrix) = cellmatrix$X
cellmatrix = cellmatrix[,-1]

metadata = read.csv("phenotype_data.csv")
rownames(metadata) = metadata$sample_id

##############get single cell and phenotype data in correct format for deconvolution

table(rownames(metadata) == colnames(cellmatrix))
scmeta = AnnotatedDataFrame(metadata)
table(rownames(scmeta@data) == colnames(cellmatrix))
scdata = ExpressionSet(assayData = as.matrix(cellmatrix), phenoData = scmeta)

##############get single cell and phenotype data in correct format for marker identification

library(Seurat)
library(limma)

cutoff = 0.0001

mcadata = CreateSeuratObject(counts = cellmatrix)
table(order(rownames(metadata)) == order(colnames(cellmatrix)))
Idents(mcadata) = metadata$absclust3
mcadata = NormalizeData(object = mcadata, verbose = FALSE)
metadata$absclust3 = as.factor(metadata$absclust3)
markerlist = list()
for (i in 1:length(levels(metadata$absclust3))){
  now = subset(FindMarkers(mcadata, ident.1 = levels(metadata$absclust3)[i], ident.2 = NULL, only.pos = TRUE), p_val_adj < cutoff)
  markerlist[[i]] = rownames(now)
  rm(now)
}
names(markerlist) = levels(metadata$absclust3)

library(reshape2)
markers = melt(markerlist)
head(markers)

for (i in 1:length(markers$value)){
  markers$count[i] = sum(markers$value == markers$value[i])}
markers = subset(markers, count == 1)
table(markers$L1)

markersnew = split(markers$value, markers$L1)

#Load my data

colhepato = genehepato.ind.fc[["counts"]]
colhepato = colhepato[order(match(rownames(colhepato), colnames(cellmatrix))),]

#run
devtools::install_github('shenorrlab/bseq-sc')

library(bseqsc)
bseqsc_config(error = TRUE)

plotCellTotals(scdata, 'absclust3','pb')

library(xbioc)
x = cpm_cell_type(scdata, clusters = 'absclust3', samples = 'pb')
ids <- intersect(unlist(markersnew), rownames(x))
x <- x[ids, , drop = FALSE]
clusters <- as.character(pVar(x, scdata@phenoData@data$absclust3))
samples <- as.character(pVar(x, scdata@phenoData@data$pb))


cpm_cell_type <- function(x, clusters, samples){
  # extract variables
  clusters <- pVar(x, clusters)
  samples <- pVar(x, samples)
  # compute total count
  tot <- colSums(exprs(x))
  # compute CPM
  cpm <- sweep(exprs(x), 2L, tot, '/')
  # re-scale with within cell type average
  sc <- as.character(paste(clusters, samples))
  tot_map <- sapply(split(tot, sc), mean)
  cpm <- sweep(cpm, 2L, tot_map[sc], '*')
  # result
  res <- x
  exprs(res) <- cpm
  res
}

library(xbioc)
x = cpm_cell_type(scdata, clusters = 'absclust3', samples = 'pb')
ids <- intersect(unlist(markersnew), rownames(x))
x <- x[ids, , drop = FALSE]
clusters <- as.character(pVar(x, scdata@phenoData@data$absclust3))
samples <- as.character(pVar(x, scdata@phenoData@data$pb))

## B <- bseqsc_basis(scdata, markersnew, clusters = 'absclust3', samples = 'pb', ct.scale = TRUE)

B = data.frame()
celltypes_keep = levels(as.factor(scdata@phenoData@data$absclust3))
genes = markers$value
genes = gsub("-","_",genes)
for (j in 1:length(celltypes_keep)){
  filter <- colnames(scdata)[scdata@phenoData@data$absclust3==celltypes_keep[j]]
  expset.filt <- scdata[,filter]
  filter <- rownames(expset.filt)[rownames(expset.filt) %in% genes]
  expset.filt <- expset.filt[filter,]
  mat_now = exprs(expset.filt)
  mat_now = mat_now[genes,]
  for (i in 1:length(rownames(mat_now))){
    B[i,j] = mean(mat_now[i,])}}
rownames(B) = genes
colnames(B) = celltypes_keep
heatmap(as.matrix(B), Colv = NA, Rowv = NA, labRow = FALSE)

fit <- bseqsc_proportions(colhepato, B, verbose = TRUE, perm=1000) #Need to add the number of permutations, otherwise it will return a p-value of 9999
props = fit$coefficients
write.csv(props, "results/single_cell_proportions.csv")

################# Make Stacked Bar Plot of Proportions
library(ggplot2)
library(reshape2)
bardf = read.csv("single_cell_proportions.csv")
bardf[1,1]="early-troph"
bardf[2,1]="mid-troph"
bardf[3,1]="late-ring"
bardf[4,1]="late-troph"
bardf[5,1]="early-shizont"
bardf[6,1]="late-shizont"
bardf[7,1]="early-ring"
bardf[8,1]="mid-shizont"

bardf$X=as.factor(bardf$X)
bardf$X = factor(bardf$X, levels=c("early-ring", "late-ring", "early-troph", "mid-troph", "late-troph", "early-shizont", "mid-shizont", "late-shizont", "F", "M"))


bardf <- melt(bardf) ##Assuming your data is a matrix. i.e. the people's names are the row and col names.
colnames(bardf) <- c("Stage","Sample_ID","Proportion")
bardf

#install.packages("wesanderson")
library(wesanderson)

ggplot(bardf, aes(x=Sample_ID, y=Proportion, fill=Stage)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Sample Name") +
  ylab("Cell Proportions") 


#################### Make correlation table of cell types
library(tidyr)

corrdf = read.csv("single_cell_proportions.csv")
corrdf[1,1]="early-troph"
corrdf[2,1]="mid-troph"
corrdf[3,1]="late-ring"
corrdf[4,1]="late-troph"
corrdf[5,1]="early-shizont"
corrdf[6,1]="late-shizont"
corrdf[7,1]="early-ring"
corrdf[8,1]="mid-shizont"

corrdf[11,1]="troph"
corrdf[12,1]="ring"
corrdf[13,1]="shizont"

for (i in 2:30){
  corrdf[11,i]=sum(corrdf[1,i],corrdf[2,i],corrdf[4,i])
  corrdf[12,i]=sum(corrdf[3,i],corrdf[7,i])
  corrdf[13,i]=sum(corrdf[5,i],corrdf[8,i],corrdf[6,i])
  }
  

corrdf <- melt(corrdf) ##Assuming your data is a matrix. i.e. the people's names are the row and col names.
colnames(corrdf) <- c("Stage","Sample_ID","Proportion")

#Make long table wide
corrdf <- spread(corrdf, Stage, Proportion)

#plot troph vs. ring

corrlm <- lm(troph ~ ring, data=corrdf)
summary(corrlm)

ggplot(data=corrdf, aes(x=ring, y=troph)) +
  geom_point() +
  geom_smooth(method="lm") +
  annotate("text", x=0.4, y=0.25, label= "R^2=0.08") +
  annotate("text", x=0.4, y=0.23, label="p-value=0.13") +
  ylab("Trophozoite") +
  xlab("Ring") 

#Plot Female vs. Ring

corrlm <- lm(F ~ ring, data=corrdf)
summary(corrlm)

ggplot(data=corrdf, aes(x=ring, y=F)) +
  geom_point() +
  geom_smooth(method="lm") +
  annotate("text", x=0.4, y=0.35, label= "R^2=0.10") +
  annotate("text", x=0.4, y=0.33, label="p-value=0.08") +
  ylab("Female") +
  xlab("Ring") 

#Plot Male vs. Female

corrlm <- lm(F ~ M, data=corrdf)
summary(corrlm)

ggplot(data=corrdf, aes(x=M, y=F)) +
  geom_point() +
  geom_smooth(method="lm") +
  annotate("text", x=0.25, y=0.35, label= "R^2=0.32") +
  annotate("text", x=0.25, y=0.33, label="p-value=0.001") +
  ylab("Female") +
  xlab("Male") 
