##### Packages #####
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggpubr)


##### PEC3.2 #####

## Load data
rna <- read.csv("Data/scRNA_2.csv", row.names = 1)
adt <- read.csv("Data/scADT_2.csv", row.names = 1)


# Create Seurat Object with RNA counts
dataset2 <- CreateSeuratObject(counts = rna, 
                              project = "PEC3")

# Add Antibody-Derived-Tags (ADTs) to the Seurat object
adt_assay <- CreateAssayObject(counts = adt)
dataset2[["ADT"]] <- adt_assay

# Mitocondrial counts already filtered in the data set
# Check for Nº features and Counts
p1 <- VlnPlot(dataset2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
p2 <- FeatureScatter(dataset2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## Normalize data: 
# RNA Normalization: ln([gene1/lib.size] * 10000 + 1)
dataset2 <- NormalizeData(dataset2)

# ADT Normalization: logarithm of the ratio between 
# the individual elements and the geometric mean of the vector
dataset2 <- NormalizeData(dataset2, normalization.method = "CLR", margin = 2, assay = "ADT")

## Find variable features
dataset2 <- FindVariableFeatures(dataset2, selection.method = "vst", nfeatures = 200)
top10 <- head(VariableFeatures(dataset2), 10)

plot1 <- VariableFeaturePlot(dataset2)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


## Scaling data (mean=0 and sd=1 for every gene across cells)
all.genes <- rownames(dataset2)
dataset2 <- ScaleData(dataset2, features = all.genes)

# Run PCA: this will run automatically using the scaled data
dataset2 <- RunPCA(dataset2)
VizDimLoadings(dataset2, dims = 1:2, reduction = "pca")
DimPlot(dataset2)

## Significant components
dataset2 <- JackStraw(dataset2, num.replicate = 100, prop.freq = 0.1)
dataset2 <- ScoreJackStraw(dataset2, dims = 1:20)
JackStrawPlot(dataset2, dims = 1:15, xmax = 1, ymax = 1)
ggsave("scRNA2 JackStraw.png", device = "png", dpi = "retina")

## Find Neighbors 
dataset2 <- FindNeighbors(dataset2, dims = 1:14)

## Cluster Cells
dataset2 <- FindClusters(dataset2, resolution = 0.3)
dataset2 <- RunUMAP(dataset2, dims = 1:14)
dataset2 <- RunTSNE(dataset2, dims = 1:14)
dataset2@meta.data$ID <- metadata$Prediction_Ind


# Find All markers
rna.markers <- FindAllMarkers(dataset2, assay = "RNA")
adt.markers <- FindAllMarkers(dataset2, assay = "ADT")

rna.markers$Pos.Neg <- ifelse(rna.markers$avg_log2FC > 0, "+", "-")
adt.markers$Pos.Neg <- ifelse(adt.markers$avg_log2FC > 0, "+", "-")
adt.markers$gene <- gsub("-AB", "", adt.markers$gene)

a <- rna.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  summarise(RNA.Markers = paste(gene, collapse = ", "))

b <- adt.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  summarise(ADT.Pos.Markers = paste(gene, collapse = ", "))


left_join(a,b) %>% 
  write.table("Markers 2.txt", 
              sep = "\t", 
              row.names = F,
              quote = F)


# Label Clusters
IDs <- dataset2@meta.data$seurat_clusters
levels(IDs) <- c("Neut",      # 0
                 "B",         # 1
                 "T",         # 2
                 "NK",        # 3
                 "EP/MkP",    # 4
                 "HSC",       # 5
                 "MP",        # 6
                 "MPP",       # 7
                 "Monocytes", # 8
                 "pDC",       # 9
                 "iB")        # 10

dataset2@meta.data$ID <- IDs %>% as.character()


## Plot UMAP, T-SNE and PCA
DimPlot(dataset2, reduction = "tsne", 
        group.by = "ID", 
        label = T, 
        label.size = 4.5,
        pt.size = 1) +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_blank(),
        legend.position = "none")

ggsave("scRNA2 tSNE.png", device = "png", dpi = "retina")

DimPlot(dataset2, reduction = "umap", 
        group.by = "ID", 
        label = T, 
        label.size = 4.5,
        pt.size = .8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(plot.title = element_blank(),
        legend.position = "none")

ggsave("scRNA2 UMAP.png", device = "png", dpi = "retina")

DimPlot(dataset2, reduction = "pca", 
        group.by = "ID", 
        label = T, 
        label.size = 4.5,
        pt.size = 1)  + 
  xlab("PCA 1") +
  ylab("PCA 2") +
  theme(plot.title = element_blank(),
        legend.position = "none")

ggsave("scRNA2 PCA.png", device = "png", dpi = "retina")


## Save Seurat Object
save(x = dataset2, file = "Dataset2.RData")




##### Pseudo-Bulk #####
## Get RNA and Protein count matrices
rna <- GetAssayData(dataset2, assay = "RNA", slot = "counts") %>% data.frame
pro <- GetAssayData(dataset2, assay = "ADT", slot = "counts") %>% data.frame

# Get metadata
metadata <- dataset2@meta.data

# Count cells per cluster
n.cell <- metadata %>% group_by(ID) %>% summarise(n = n()) %>% data.frame %>% filter(., complete.cases(.))

# Keep Clusters with over n cells
n <- 400
kp <- n.cell$ID[n.cell$n > n]

# Filter Metadata File with selected cells
metadata.filt <- filter(metadata, ID %in% kp)


# Merge Data sets
# Change Row Names on Protein data to ENSEMBL Gene ID
new.rows <- read.table("PEC3/AB.info.txt", sep = "\t", header = T)
pro <- pro[new.rows$Antibody,]
rownames(pro) <- new.rows$Gene.stable.ID
pro <- pro[complete.cases(pro),]

# Change Row Names on RNA data to ENSEMBL Gene ID
new.rows <- read.table("PEC3/Name_2_ID.txt", sep = "\t", header = T)
rna <- rna[new.rows$Gene.name,]
rownames(rna) <- new.rows$Gene.stable.ID
rna <- rna[complete.cases(rna),]

rna.f <- filter(data.frame(rna), rownames(rna) %in% rownames(pro))
pro.f <- filter(data.frame(pro), rownames(pro) %in% rownames(rna))

rna <- rna.f[rownames(pro.f),] %>% data.frame()
pro <- pro.f %>% data.frame()


## Compute rna.bulk counts (No down-sampling see test below)
rna.bulk <- tapply(rownames(metadata.filt), metadata.filt$ID, function(a) {
  cells <- rna[,colnames(rna) %in% a]
  rowSums(cells)
})

pro.bulk <- tapply(rownames(metadata.filt), metadata.filt$ID, function(a) {
  cells <- pro[,colnames(pro) %in% a]
  rowSums(cells)
})

# Transform into matrix
rna.bulk <- do.call(cbind, rna.bulk)
pro.bulk <- do.call(cbind, pro.bulk)

## Normalize data
rna2 <- edgeR::cpm(rna.bulk) %>% log1p()
pro2 <- edgeR::cpm(pro.bulk) %>% log1p()

## Scale to equal variance
rna2 <- scale(rna2) + mean(colMeans(rna1))
pro2 <- scale(pro2) + mean(colMeans(pro1))



## Compute CF
cf1 <- pro1/rna1
cf2 <- pro2/rna2

# Compute median CF
cf1 <- apply(cf1, 1, function(a) median(a, na.rm = T))
cf2 <- apply(cf2, 1, function(a) median(a, na.rm = T))


# Correct RNA
rna2.corr1 <- sweep(rna2, 1, cf1[names(cf2)], "*")
rna2.corr2 <- sweep(rna2, 1, cf2, "*")


# Apply CF
## Test Correlations
corrs1 <- c()
corrs2 <- c()
corrs3 <- c()



for(i in seq_along(colnames(rna2))){
  # Get RNA and Protein values from tissue i
  t1.rna <- rna2[,i]
  t1.pro <- pro2[,i]
  t1.rna.corr1 <- rna2.corr1[,i]
  t1.rna.corr2 <- rna2.corr2[,i]
  
  # Compute Correlations
  kp <- !is.na(t1.rna.corr1)
  
  corrs1[i] <- cor(t1.rna[kp], t1.pro[kp], method = "s")
  corrs2[i] <- cor(t1.rna.corr1[kp], t1.pro[kp], method = "s")
  corrs3[i] <- cor(t1.rna.corr2[kp], t1.pro[kp], method = "s")
}

## Store correlations before and after correction.
names(corrs1) <- names(corrs2) <- names(corrs3) <- colnames(rna2)
corrs.sc2 <- cbind("Uncorrected" = corrs1, "Corrected.CF1" = corrs2, "Corrected.CF2" = corrs3)


## Boxplot
corrs.plot <- gather(data.frame(corrs.sc2))

corrs.plot$key <- factor(corrs.plot$key, 
                         levels = c("Uncorrected", "Corrected.CF1", "Corrected.CF2"),
                         labels = c("Uncorrected", "Corrected CF1", "Corrected CF2"),
                         ordered = T)

ggplot(corrs.plot, aes(x = key, y = value)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = key),
              position = position_jitter(.2), 
              size = 2,
              alpha = .7) +
  ylab("Rho") +
  ggtitle("Comparison") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 15),
        legend.position = "null",
        plot.title = element_text(hjust = .5)) +
  scale_color_manual(values = c("#7F7F7F", "#5B9BD5", "#5B9BD5"))


ggsave(filename = paste0("Single Cell Scatter Plots/", "SC2 Boxplot",".png"), dpi = "retina")


## Scaterplots
# Before Correction
for (i in 1:9) {
  ggplot(NULL, aes(x = rna2[,i], y = pro2[,i])) +
    theme_test() +
    geom_point(colour = "#7F7F7F") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("ADT (Predicted)") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna2)[i], " | ","rho =", round(corrs.sc2[i,1], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Single Cell Scatter Plots/", "Uncorrected2.scCF2.", i,".png"), dpi = "retina")
}

# After Correction
for (i in 1:9) {
  # CF1
  ggplot(NULL, aes(x = rna2.corr1[,i], y = pro2[,i])) +
    theme_test() +
    geom_point(colour = "#5B9BD5") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("ADT (Predicted)") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna2)[i], " | ","rho =", round(corrs.sc2[i,2], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Single Cell Scatter Plots/", "Corrected.scCF1.", i,".png"), dpi = "retina")
  
  
  # CF2
  ggplot(NULL, aes(x = rna2.corr2[,i], y = pro2[,i])) +
    theme_test() +
    geom_point(colour = "#5B9BD5") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("ADT (Predicted)") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna2)[i], " | ","rho =", round(corrs.sc2[i,3], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Single Cell Scatter Plots/", "Corrected.scCF2.", i,".png"), dpi = "retina")
}

