############## Immune cell composition PCA ############
library(factoextra)

setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR/")
immune = read.table("Simons_etal_2019_immune_components.txt",
                    header=TRUE)
#immune$Sample[immune$Sample == "RC106"] = "RC106R"
immune = immune[order(immune$Sample),] # Put in order
rownames(immune) = immune$Sample
immune.ind = immune[,-1:-2]

# pca by individual
pca_res <- prcomp(immune.ind) 

#pca by cell type
immune.cell <- data.frame(t(immune.ind))
pca_cell <- prcomp(immune.cell)

# plot them 
plot.pca.ind = data.frame(pca_res$x) # PCs vs. Individuals
plot(PC2 ~ PC1, data = plot.pca.ind)

plot.pca.cell = data.frame(pca_cell$x) # PCs vs. Cell type
plot(PC2 ~ PC1, data = plot.pca.cell)

#parasitemia_sub = subset(parasitemia, Individual_ID != "RTi13" | Time_point != 2)
#parasitemia_sub = subset(parasitemia_sub, Individual_ID != "RUn13" | Time_point != 2)
library(factoextra)
# Look at eigenvalues
fviz_eig(pca_res) # by individual
fviz_eig(pca_cell) # by cell
get_eigenvalue(pca_cell)

get_pca_ind(pca_cell)
#Graph of individuals
fviz_pca_ind(pca_res,
             col.ind = "contrib", # Color by the contribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) +
  labs(x = "PC1", y = "PC2")

#by quality of cells
fviz_pca_ind(pca_cell,
             col.ind = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) +
  labs(x = "PC1", y = "PC2", title = "Cells - PCA")


# look at which PC parasitemia contributes to
fviz_contrib(pca_res, choice = "ind", axes = 1)
fviz_contrib(pca_cell, choice = "ind", axes = 1)

############### Results for Variables ############
res.var <- get_pca_var(pca_cell)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Extract and save results
write.table(res.var$contrib, file="edgeR_results/PCA_results",
            sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
