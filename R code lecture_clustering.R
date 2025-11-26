
# Find path to your files by: Session-> Set Working Directory -> Choose Directory
# Press open, and copy and paste the path printed in your console

setwd("Your path to the files")
#setwd("C:/NGS/scRNAseq") # Example

#Create data in R for clustering
rm(list=ls()) 
A=c(70,130 ,100 ,180, 190)
B=c(86,147 ,90 ,170,215)
C=c(78, 137 ,100, 181, 192)
D=c(60, 121 ,90 ,171 ,215)
E=c(75,134 ,106 , 185, 200)
df=data.frame(A,B,C,D,E)
variables=c("DBP","SBP","Weight","Height","Cholesterol")
rownames(df)=variables
df

#Calculate distance matrix in R
h=dist(df, method = "euclidean")
h

#Create dendrogram in R
hp=hclust(h,method="average")
plot(hp)

#Identify clusters
cutree(hp, h=150)
cutree(hp, h=100)

# Linkage function – mean (UPGMA)
hp=hclust(h,method="average")
plot(hp)

# Single linkage function
hp=hclust(h,method="single")
plot(hp)

#Complete linkage function
hp=hclust(h,method="complete")
plot(hp)

# Example of a bootstrap
v=1:10 
sample(v,replace=TRUE)

#Select a linkage function that is most robust
v=1:5
s1=sample(v,replace=TRUE)
boot=df[,s1]
boot
h=dist(boot, method ="euclidean")
hp=hclust(h,method="average")
plot(hp)

# Use the pvclust
library(pvclust)
set.seed(1970) #For reproducibility
fit = pvclust(t(df), method.hclust="average",     
              method.dist="euclidean",nboot=100)
mean(fit$edges$bp)
plot(fit)

#Standardize the variables
sdf=t(scale(t(df)))
sdf

# Cluster based on standardized variables
h=dist(sdf, method = "euclidean")
hp=hclust(h,method="average")
plot(hp)

# Correlation matrix
round(cor(t(df)),2) 

# Compute distance matrix based on correlations
cor_dist=1-cor(t(df)) 
h=as.dist(cor_dist) 
round(h,2) 

# Cluster based on correlations
hp=hclust(h,method="average") 
plot(hp,ylab="1-Pearson's correlation") 

#Transpose the data
tdf=t(df) 

# Cluster based on transposed data
h=dist(tdf, method = "euclidean")
hp=hclust(h,method="average")
plot(hp) 

# Heatmaps
library(pheatmap)
pheatmap(as.matrix(df),scale="none",clustering_method="average")

# Scale the rows
library(pheatmap)
pheatmap(as.matrix(df),scale="row",clustering_method="average")



# Load the drug data
rm(list=ls())
df=read.csv("Cytokines.csv",header=TRUE,sep=";",dec=",")
df[1:4,1:4]

# Extract only dugs (not the controls) and log the data
library(dplyr)
library(tibble)
Drugs = df %>%
  filter(!Treatment %in% c("ConA_DMSO", "ConA")) %>%
  column_to_rownames(var = "Treatment")

log10.Drugs=t(log10(Drugs))
log10.Drugs[1:4,1:3]

# Correlation matrix
round(cor(log10(Drugs)),2)

# Clustering
library(pvclust)
fit = pvclust(log10(Drugs),method.hclust="average",
              method.dist="cor",nboot=1000)
plot(fit)

# Heatmap
library(pheatmap) 
pheatmap(as.matrix(t(log10.Drugs)),
         scale="column", 
         clustering_distance_cols="correlation", 
         clustering_method = "average", 
         fontsize_row=3)



#########################
# Kmeans clustering
rm(list=ls()) 
x=c(1,2,4,5,7,8,10,11,9,12)
y=c(9,7,3,3,2,4,11,9,10,8)
df=data.frame(x,y)
plot(x,y,cex=2,pch=16)
text(x,y,1:10,cex=0.8,col="white")

k_clust = kmeans(df, 3,nstart = 25) 
k_clust$cluster

# Plot
plot(x,y,cex=2,pch=16,col=k_clust$cluster)
text(x,y,1:10,cex=0.8,col="white")

# The elbow method
n=8
WCSS=NULL
for (i in 1:n) {
  WCSS[i]=kmeans(df, i, nstart = 25)$tot.withinss
}

plot(1:n, WCSS,type="b", pch = 19, 
     xlab="Number of clusters (k)",
     ylab="Total within cluster sum of squares")

#################################
# Graph based clustering (Louvain) 
#################################

## kNN graph
library(igraph)   # for graph building
x=c(0.08,0.22,0.12,0.44,0.65,0.70)
y=c(0.8,0.78,0.5,0.35,0.1,-0.3)
df <- data.frame(x,y)
(dist_mat <- as.matrix(dist(df)))

k <- 2 # number of neighbors = 2 (+1 for the main diagonal) 
n <- nrow(dist_mat)
adj_matrix <- matrix(0, n, n)

for (i in 1:n) {
  # Including self (distance = 0)
  neighbors <- order(dist_mat[i, ])[1:(k+1)]
  adj_matrix[i, neighbors] = 1
}
diag(adj_matrix)=0 # Set self to 0
adj_matrix

library(igraph)   # for graph building
g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed")
plot(g, edge.arrow.size = 0.6,vertex.size = 20,vertex.label.cex = 0.9)

# Weighted SNN graph (Seurat)
library(Seurat) # Version 5.3.0
#Seurat (includes self node as neighbor)
obj <- FindNeighbors(dist(df) ,k.param = k+1,prune.SNN=0.25)
SNNS <- as.matrix(obj$snn)
diag(SNNS)=0 # Skip self node as neighbor
SNNS # print weights
library(igraph)
g = graph_from_adjacency_matrix(SNNS,mode = "undirected", weighted=TRUE)
plot(g)
cluster_louvain(g)
##################################################################
## Optional Seurat clustering
##################################################################

library(Seurat)
rm(list=ls())
set.seed(42)
scRNAseq_data=read.csv("scRNAseq_data.csv", row.names = 1)
pbmc = CreateSeuratObject(counts = scRNAseq_data,project = "my_project",min.cells = 0,min.features = 0)
pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc= FindVariableFeatures(pbmc, selection.method = "dispersion", nfeatures = 2000) 
pbmc= ScaleData(pbmc) 
pbmc= RunPCA(object = pbmc, features = VariableFeatures(object = pbmc)) 
pbmc= FindNeighbors(object = pbmc, dims = 1:10) # Construct SNN graph
pbmc= FindClusters(pbmc, resolution = 0.5) # Set resolution to 0.5
snn_graph = pbmc@graphs$RNA_snn # Extract SNN graph
g = graph_from_adjacency_matrix(as.matrix(snn_graph), mode = "undirected", weighted = TRUE, diag = FALSE)
V(g)$cluster = as.factor(pbmc$seurat_clusters)
V(g)$color = V(g)$cluster  # igraph will color by factor
plot(g,layout = layout_with_fr, vertex.size = 4,
     vertex.label = NA,
     edge.width = 0.1,  
     vertex.color = V(g)$color)
pbmc = RunUMAP(object = pbmc, dims = 1:10) # Perform UMAP
DimPlot(pbmc, label = TRUE, reduction = "umap") # Plot UMAP

##################################################################
## Optional scRNAseq Louvain (SNN) clustering with manual code
##################################################################
rm(list=ls())
# 1 Read in scRNAseq data
scRNAseq_counts = read.csv("scRNAseq_data.csv", row.names = 1, check.names = FALSE)
ncells=ncol(scRNAseq_counts) # Number of cells
# 2 Normalize
norm_counts=scRNAseq_counts # Copy matrix
for (j in 1:ncells) {
  norm_counts[, j] = (scRNAseq_counts[, j] / sum(scRNAseq_counts[, j])) * 10000
}
log_norm_counts=log(norm_counts+1)
# 3 Identify HVG
gene_vars = apply(norm_counts, 1, var) # Calculate variance of genes
gene_means = apply(norm_counts, 1, mean) # Calculate mean of genes
Dispersion=gene_vars/gene_means
top2000names = names(sort(Dispersion, decreasing = TRUE))[1:2000]
log_norm_counts_top2000=log_norm_counts[top2000names,] # Extract top 2000 genes
# 4 Scale genes
scaled_data_top2000 = t(scale(t(log_norm_counts_top2000)))
# 5 PCA
pca_res = prcomp(t(scaled_data_top2000), center = F, scale. = F)
pca_embed = pca_res$x[, 1:10]
# 6. Louvain SNN clustrering
library(bluster)
library(igraph)
# Seurat k = 20 -> 19
SNN <- makeSNNGraph(pca_embed, k = 19,type="jaccard")# jaccard/number
A <- as_adjacency_matrix(SNN, attr = "weight", sparse = FALSE)
# Prune
A[A<1/15]=0
g <- graph_from_adjacency_matrix(A, mode = "undirected",weighted=T)
louvain=cluster_louvain(g, resolution = 0.8)
cluster=factor(membership(louvain)) 
# Plot SNN graph
# Create a plot that shows the clusters in the SNN graph
V(g)$cluster = cluster
V(g)$color = cluster
plot(g,layout = layout_with_fr, vertex.size = 4,
     vertex.label = NA,
     edge.width = 0.1,  
     vertex.color = V(g)$color)


# 7. UMAP
library(umap)
set.seed(43)
umap_res <- umap(pca_embed, 
                 n_neighbors = 30, min_dist=0.3)
UM=data.frame(Dim1 = umap_res$layout[, 1],
              Dim2 = umap_res$layout[, 2],
              Cluster = cluster)

ggplot(UM, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.6) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "UMAP single cell data",color="Cluster",
       x = "UMAP Dimension 1", y = "UMAP Dimension 2")



























