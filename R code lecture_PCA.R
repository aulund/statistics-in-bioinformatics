
A=matrix(c(67,77,89,170,173,179,120,123,130),3,3)
A
A[3,2]
A[3,1:3] # Print all values in row 3
A[3,] # Print all values in row 3 


A>100 # Check which elements that have a value > 100


A[,1:2] # print all rows in columns 1 and 2
A[,c(1,3)] # print all rows in columns 1 and 3


A
t(A) # Transpose


# Vectors
V=matrix(c(2,3,5,1),4,1)
V # Print
t(V)

# Matrix addition
A = matrix(c(1,3,2,4), 2,2)
B = matrix(2, 2,2)
A + B

#scalar multiplication.
3*A

# Element-wise multiplication in R 
A = matrix(c(1,3,2,4), 2,2)
B = matrix(c(5,7,6,8), 2,2)
A * B

# Matrix multiplication (dot product) in R
A = matrix(c(1,3,2,4), 2,2)
B = matrix(c(5,7,6,8), 2,2)
A %*% B

# Eigenvectors and Eigenvalues in R
A=matrix(c(1,1,2,0),2,2)
eigen(A)

# PCA
# 1. Center the data
SBP=c(126,128,128,130,130,132)
DBP=c(78,80,82,82,84,86)
data=data.frame(SBP,DBP)
cdata=scale(data,center=T,scale=F)
cdata # print

# 2. Calculate the covariance matrix
COV=cov(cdata)
COV # print

# 3. Calculate the eigenvectors of the covariance matrix 
Eig=eigen(COV)
Eig # print

# 4. Order the eigenvectors
# The eginvectors are already ordered

# 5. Calculate the principal components
PCs=cdata%*%Eig$vectors
round(PCs,2) # print

# Plot the PCs
plot(PCs[,1],PCs[,2],ylim=c(-5,5),pch=16, cex=3,col="darkgreen")
abline(h=0)
abline(v=0)
# The prcomp function
pca=prcomp(data,center=T,scale=F)
round(pca$x,2)

# Compute the variances of the PCs
pca$sdev^2

# Compute proportion of variance explained
summary(pca)

# Download the 3 data files to the same folder
# Set directory in RStudio: Session -> Set Working Directory -> Choose Directory...

# Make sure you are in the correct folder by listing the files:
list.files()

# Read in the bulk RNAseq data in R
Bulk_RNAseq_counts <- read.csv("/Users/aulund/Documents/GitHub/statistics in bioinf/CVID_bioinfo.csv", header=TRUE, sep=";", row.names = 1)
head(Bulk_RNAseq_counts)

Group=factor(rep(c("CVID","HD"),c(9,9))) # Define groups

t_Bulk_RNAseq_counts=t(Bulk_RNAseq_counts) # Transpose
t_Bulk_RNAseq_counts[1:4,1:4] # print

# PCA
pca_res = prcomp(t_Bulk_RNAseq_counts, center = TRUE, scale. = TRUE)


# Plot PC1 vs PC2
PCA = pca_res$x[, 1:2]
plot(PCA,pch=16,cex=2,col=Group,main="PCA")
text(PCA[,1],PCA[,2],1:18,cex=0.6,col="white")

# Bar chart of first 10 components
var_explained = (pca_res$sdev)^2
prop = var_explained / sum(var_explained)
barplot(prop[1:10],xlab = "PCs",
        ylab = "%Variance",col = "skyblue",names.arg = 1:10)


# Plot with ggplot2
library(ggplot2)
ggplot(PCA, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_text(aes(PC1, PC2,label=1:length(PC1)),size=4,color="black") +
  theme_minimal() +
  labs(x = paste("PC1 (",round(prop[1]*100) ,"%)", sep=""),
       y = paste("PC2 (",round(prop[2]*100) ,"%)", sep=""),
       title = "Bulk RNAseq 25 genes")

# We will now try PCA on scRNAseq data
# 1. Read in the sc data
rm(list=ls())
scRNAseq_counts = read.csv("/Users/aulund/Documents/GitHub/statistics in bioinf/scRNAseq_data.csv", row.names = 1, check.names = FALSE)
dim(scRNAseq_counts)
scRNAseq_counts[1:12,1:3]

# 2. Normalize
ncells=ncol(scRNAseq_counts) # Number of cells
norm_counts=scRNAseq_counts # Copy matrix
for (j in 1:ncells) {
  norm_counts[, j] = (scRNAseq_counts[, j] / sum(scRNAseq_counts[, j])) * 10000
}
log_norm_counts=log(norm_counts+1)
log_norm_counts[1:15,1:2] # Print

# Variance vs mean on normalized counts across the cells
gene_vars = apply(norm_counts, 1, var) 
gene_means = apply(norm_counts, 1, mean)
plot(gene_means,gene_vars,log="xy")

# Find HVG (the ones with large dispersion)
Dispersion=gene_vars/gene_means
top2000names = names(sort(Dispersion, decreasing = TRUE))[1:2000]
points(gene_means[top2000names],gene_vars[top2000names],cex=0.5,col="red")
log_norm_counts_top2000=log_norm_counts[top2000names,] # Extract top 2000 genes
dim(log_norm_counts_top2000)


# 4. Scaling (standardize values of HVG)
scaled_data_top2000 = t(scale(t(log_norm_counts_top2000)))
scaled_data_top2000[1:4,1:2]


# 5. PCA (this may take a minute)
pca_res = prcomp(t(scaled_data_top2000), center = F, scale. = F)

# Plot variance explained
var_explained = (pca_res$sdev)^2
var_explained_ratio = var_explained / sum(var_explained)
# Bar chart of first 100 components
barplot(var_explained_ratio[1:30],xlab = "PCs",
        ylab = "Proportion of Variance",
        col = "skyblue",names.arg = 1:30)

# Extract the 10 first components
pca_embed = pca_res$x[, 1:10]

# 6. Clustering with K-means
# Identify 7 clusters based on the first 10 PCs
km = kmeans(pca_embed, centers = 7, nstart = 200)
cluster=factor(km$cluster)
cluster[1:3]


# Plot cell clusters with PCA
# PC1 vs PC2
plot(pca_embed[,1:2],pch=16,cex=0.5,col=cluster,main="PCA")
# PC2 vs PC3
plot(pca_embed[,c(2,3)],pch=16,cex=0.5,col=cluster,main="PCA")

# t-SNE
library(Rtsne)
set.seed(9)
tsne_res <- Rtsne(pca_embed, dims = 2)
plot(tsne_res$Y,pch=16,cex=0.6,
     col=cluster,main="t-SNE")

# UMAP
library(umap)
set.seed(43)
umap_res <- umap(pca_embed, 
                 n_neighbors = 30, min_dist=0.3)
plot(umap_res$layout,pch=16,
     cex=0.3,col=cluster,main="UMAP")

# UMAP plot in ggplot2
# Nicer plot with ggplot2
UM=data.frame(Dim1 = umap_res$layout[, 1],
              Dim2 = umap_res$layout[, 2],
              Cluster = cluster)

ggplot(UM, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.6) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "UMAP single cell data",color="Cluster",
       x = "UMAP Dimension 1", y = "UMAP Dimension 2")

# Identify cells (annotate cells)
# Choose a gene: CD3D = T cells, MS4A1 = B cells, CD14 = Monocytes, , KLRD1 = NK-cells
gene = "MS4A1"
expr = as.numeric(log_norm_counts[gene, ])
cols = colorRampPalette(c("lightgrey", "blue"))(100)[cut(expr, breaks = 100)]
plot(umap_res$layout,pch=16,cex=0.5,col=cols,main="UMAP")



# Drug data in R
rm(list=ls())
df=read.csv("/Users/aulund/Documents/GitHub/statistics in bioinf/Cytokines.csv",header=TRUE,sep=";",dec=",")
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

# PCA
pca=prcomp(log10.Drugs,scale=TRUE)
# Bar chart of first 10 components
var_explained <- (pca$sdev)^2
var_explained_ratio <- var_explained / sum(var_explained)
barplot(var_explained_ratio[1:10],
        xlab = "Principal Component",
        ylab = "Proportion of Variance",
        col = "skyblue",names.arg = 1:10)

# Plot PC1 vs PC2
PCA = pca$x[, c(1,2)]
plot(PCA,pch=16,cex=2,main="PCA",col="lightblue")
text(PCA[,1],PCA[,2],rownames(PCA),cex=0.6)

# Plot PC1 vs PC3
PCA = pca$x[, c(1,3)]
plot(PCA,pch=16,cex=2,main="PCA",col="lightblue")
text(PCA[,1],PCA[,2],rownames(PCA),cex=0.6)

