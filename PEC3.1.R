##### Packages #####
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggpubr)


##### PEC3.1 #####

## Load data
rna <- read.csv("Data/scRNA_1.csv", row.names = 1)
adt <- read.csv("Data/scADT_1.csv", row.names = 1)

# Create Seurat Object with RNA counts
dataset <- CreateSeuratObject(counts = rna, project = "PEC3")

# Add Antibody-Derived-Tags (ADTs) to the Seurat object
adt_assay <- CreateAssayObject(counts = adt)
dataset[["ADT"]] <- adt_assay

# Mitocondrial counts already filtered in the dataset
# Check for Nº features and Counts
p1 <- VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
p2 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## Normalize data: 
# RNA Normalization: ln([gene1/lib.size] * 10000 + 1)
dataset <- NormalizeData(dataset)

# ADT Normalization: logarithm of the ratio between 
# the individual elements and the geometric mean of the vector
dataset <- NormalizeData(dataset, normalization.method = "CLR", margin = 2, assay = "ADT")

## Find variable features
dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 200)
top10 <- head(VariableFeatures(dataset), 10)
plot1 <- VariableFeaturePlot(dataset)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


## Scaling data (mean=0 and sd=1 for every gene across cells)
all.genes <- rownames(dataset)
dataset <- ScaleData(dataset, features = all.genes)

# Run PCA: this will run automatically using the scaled data
dataset <- RunPCA(dataset)
VizDimLoadings(dataset, dims = 1:2, reduction = "pca")
DimPlot(dataset)

## Significant components
dataset <- JackStraw(dataset, num.replicate = 100, prop.freq = 0.1)
dataset <- ScoreJackStraw(dataset, dims = 1:20)
JackStrawPlot(dataset, dims = 1:15, xmax = 1, ymax = 1)
ggsave("scRNA1 JackStraw.png", device = "png", dpi = "retina")


## Run embedings (UMAP and TSNE)
dataset <- RunUMAP(dataset, dims = 1:14)
dataset <- RunTSNE(dataset, dims = 1:14)

## Find Neighbors 
dataset <- FindNeighbors(dataset, dims = 1:14)

## Cluster Cells
dataset <- FindClusters(dataset, resolution = 0.4)

## Find All markers
rna.markers <- FindAllMarkers(dataset, assay = "RNA")
adt.markers <- FindAllMarkers(dataset, assay = "ADT")

a <- rna.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  summarise(RNA.Markers = paste(gene, collapse = ", "))

b <- adt.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  summarise(ADT.Pos.Markers = paste(gene, collapse = ", "))


left_join(a,b) %>% 
  write.table("Markers 1.txt", 
              sep = "\t", 
              row.names = F,
              quote = F)


## Relabel Clusters

IDs <- dataset@meta.data$seurat_clusters
levels(IDs) <- c("Neut",      # 0
                 "NK",        # 1
                 "B",         # 2
                 "Mono",      # 3
                 "T-CD8",     # 4
                 "T-CD4",     # 5
                 "EP/MkP",    # 6
                 "HSC/MPP",   # 7
                 "MP",        # 8
                 "CD1c DC",   # 9
                 "pDC",       # 10
                 "Myelo")     # 11

dataset@meta.data$ID <- IDs %>% as.character()


## Plot UMAP, T-SNE and PCA
DimPlot(dataset, reduction = "tsne", 
        group.by = "ID", 
        label = T, 
        label.size = 4.5,
        pt.size = 1) +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  theme(plot.title = element_blank(),
        legend.position = "none")

ggsave("scRNA1 tSNE.png", device = "png", dpi = "retina")

DimPlot(dataset, reduction = "umap", 
        group.by = "ID", 
        label = T, 
        label.size = 4.5,
        pt.size = .8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(plot.title = element_blank(),
        legend.position = "none")

ggsave("scRNA1 UMAP.png", device = "png", dpi = "retina")

DimPlot(dataset, reduction = "pca", 
        group.by = "ID", 
        label = T, 
        label.size = 4.5,
        pt.size = 1)  + 
  xlab("PCA 1") +
  ylab("PCA 2") +
  theme(plot.title = element_blank(),
        legend.position = "none")

ggsave("scRNA1 PCA.png", device = "png", dpi = "retina")



## Save Seurat Object
save(x = dataset, file = "Dataset.RData")





##### Pseudo-Bulk #####
## Get RNA and Protein count matrices
rna <- GetAssayData(dataset, assay = "RNA", slot = "counts") %>% data.frame
pro <- GetAssayData(dataset, assay = "ADT", slot = "counts") %>% data.frame

## Get Metadata
metadata <- dataset@meta.data

# Count cells per cluster
n.cell <- metadata %>% group_by(ID) %>% summarise(n = n()) %>% data.frame

# Keep Clusters with over n cells
n <- 400
kp <- n.cell$ID[n.cell$n > n]

# Filter Metadata File with selected cells
metadata.filt <- filter(metadata, ID %in% kp)

# Merge Datasets
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

# Select genes present in both, rna and pro, datasets.
rna.f <- filter(data.frame(rna), rownames(rna) %in% rownames(pro))
pro.f <- filter(data.frame(pro), rownames(pro) %in% rownames(rna))

rna <- rna.f[rownames(pro.f),] %>% data.frame()
pro <- pro.f %>% data.frame()

## Compute rna.bulk counts
rna.bulk <- tapply(rownames(metadata.filt), metadata.filt$ID, function(a) {
  cells <- dplyr::select(rna, all_of(a))
  rowSums(cells)
})


pro.bulk <- tapply(rownames(metadata.filt), metadata.filt$ID, function(a) {
  cells <- dplyr::select(pro, all_of(a))
  rowSums(cells)
})

# Transform into matrix
rna.bulk <- do.call(cbind, rna.bulk)
pro.bulk <- do.call(cbind, pro.bulk)

## Normalize Gene counts with CPM from edgeR
rna1 <- edgeR::cpm(rna.bulk) %>% log1p
pro1 <- edgeR::cpm(pro.bulk) %>% log1p

rna1 <- scale(rna1) + mean(colMeans(rna1))
pro1 <- scale(pro1) + mean(colMeans(pro1))

## Scale to match previous dataset
rna.m <- mean(colMeans(rna.pro.filt)[01:29])
pro.m <- mean(colMeans(rna.pro.filt)[30:58])

rna1 <- scale(rna1) + rna.m
pro1 <- scale(pro1) + pro.m

## Test Correlations
cf.sc <- rowMeans(cf, na.rm = T)
cf.sc <- cf.sc[rownames(rna1)]

corrs1 <- c()
corrs2 <- c()
rna.corr <- list()

for(i in seq_along(colnames(rna1))){
  # Get RNA and Protein values from tissue i
  t1.rna <- rna1[,i] %>% unlist
  t1.pro <- pro1[,i] %>% unlist
  
  # Remove 0 values to avoid correlation artifacts
  corrs1[i] <- cor(t1.rna, t1.pro, method = "s")
  
  # Predict Protein values from RNA levels using the
  # information from the other 28 tissues.
  t1.rna.corr <- t1.rna * cf.sc
  
  # Store values
  rna.corr[[i]] <- t1.rna.corr
  
  # Remove NA values to avoid correlation artifacts
  kp <- !is.na(cf.sc)
  corrs2[i] <- cor(t1.rna.corr[kp], t1.pro[kp], method = "s")
}

## Create data set with Protein predicted data
rna1.corr <- do.call(cbind, rna.corr)

## Store correlations before and after correction.
names(corrs1) <- names(corrs2) <- colnames(rna1)
corrs <- cbind("Uncorrected" = corrs1, "Corrected" = corrs2)

# Compare correlations before and after:
t.test(corrs1, corrs2, paired = T)

## Boxplot
corrs.plot <- gather(data.frame(corrs))
corrs.plot$key <- factor(corrs.plot$key, 
                         levels = c("Uncorrected", "Corrected"), 
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
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "null",
        plot.title = element_text(hjust = .5)) +
  scale_color_manual(values = c("#7F7F7F", "#5B9BD5"))

ggsave(filename = paste0("Single Cell Scatter Plots/", "Bulk Boxplot",".png"), dpi = "retina")

## Scaterplots
# Remove NAs from plots
kp <- !is.na(cf.sc)

for (i in 1:8) {
  ggplot(NULL, aes(x = rna1[kp,i], y = pro1[kp,i])) +
    theme_test() +
    geom_point(colour = "#7F7F7F") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("RNA") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna1)[i], " | ","rho =", round(corrs[i,1], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
   ggsave(filename = paste0("Single Cell Scatter Plots/", "Uncorrected.", i,".png"), dpi = "retina")
}


# After Correction
for (i in 1:8) {
  ggplot(NULL, aes(x = rna1.corr[kp,i], y = pro1[kp,i])) +
    theme_test() +
    geom_point(colour = "#5B9BD5") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("ADT (Predicted)") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna1)[i], " | ","rho =", round(corrs[i,2], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Single Cell Scatter Plots/", "Corrected.BulkCF.", i,".png"), dpi = "retina")
}






##### Create CF from single cell Data #####
# Select proportion of cells selected for training data
p <- 0.7

# Compute bulk counts
bulk.data <- tapply(rownames(metadata.filt), metadata.filt$ID, function(a) {
  
  # Select a random number of cells for the training set
  # Since random numbers are involved we use a seed generator to
  # ensure reproducibility. 
  set.seed(123)
  n <- length(a)
  train <- sample(1:n, round(n*p))
  
  # Select cells from the same cluster
  cells.rna <- select(rna, all_of(a))
  cells.pro <- select(pro, all_of(a))
  
  # Extract train
  rna.train <- rowSums(cells.rna[,train])
  pro.train <- rowSums(cells.pro[,train])
  
  # Extract test
  rna.test <- rowSums(cells.rna[,-train])
  pro.test <- rowSums(cells.pro[,-train])
  
  # Return a list with the values
  list(rna.train, rna.test,
       pro.train, pro.test)
  }
)


## Normalize data
rna.train.norm <- sapply(bulk.data, function(a) a[[1]])  %>% edgeR::cpm() %>% log1p()
pro.train.norm <- sapply(bulk.data, function(a) a[[3]]) %>% edgeR::cpm() %>% log1p()

rna.test.norm <- sapply(bulk.data, function(a) a[[2]]) %>% edgeR::cpm() %>% log1p()
pro.test.norm <- sapply(bulk.data, function(a) a[[4]]) %>% edgeR::cpm() %>% log1p()


rna1.train <- scale(rna.train.norm) + mean(colMeans(rna.train.norm))
pro1.train <- scale(pro.train.norm) + mean(colMeans(pro.train.norm))

rna1.test <- scale(rna.test.norm) + mean(colMeans(rna.test.norm))
pro1.test <- scale(pro.test.norm) + mean(colMeans(pro.test.norm))


## Compute CF
cf <- pro1.train/rna1.train

# Compute median CF
cf <- apply(cf, 1, function(a) median(a, na.rm = T))

# Correct RNA
rna1.test.corr <- sweep(rna1.test, 1, cf, "*")

# Apply CF
## Test Correlations
corrs1 <- c()
corrs2 <- c()

for(i in seq_along(colnames(rna1))){
  # Get RNA and Protein values from tissue i
  t1.rna <- rna1.test[,i]
  t1.pro <- pro1.test[,i]
  t1.rna.corr <- rna1.test.corr[,i]
  
  # Compute Correlations
  corrs1[i] <- cor(t1.rna, t1.pro, method = "s")
  corrs2[i] <- cor(t1.rna.corr, t1.pro, method = "s")
}

## Store correlations before and after correction.
names(corrs1) <- names(corrs2) <- colnames(rna1)
corrs.sc <- cbind("Uncorrected" = corrs1, "Corrected" = corrs2)

# Compare correlations before and after:
t.test(corrs1, corrs2, paired = T)

## Boxplot
corrs.plot <- gather(data.frame(corrs.sc))

corrs.plot$key <- factor(corrs.plot$key, 
                         levels = c("Uncorrected", "Corrected"),
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
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "null",
        plot.title = element_text(hjust = .5)) +
  scale_color_manual(values = c("#7F7F7F", "#5B9BD5"))

ggsave(filename = paste0("Single Cell Scatter Plots/", "SC Boxplot",".png"), dpi = "retina")


## Scaterplots
# Before Correction
for (i in 1:8) {
  ggplot(NULL, aes(x = rna1.test[,i], y = pro1.test[,i])) +
    theme_test() +
    geom_point(colour = "#7F7F7F") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("ADT (Predicted)") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna1)[i], " | ","rho =", round(corrs.sc[i,1], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Single Cell Scatter Plots/", "Uncorrected.scCF.", i,".png"), dpi = "retina")
}

# After Correction
for (i in 1:8) {
  ggplot(NULL, aes(x = rna1.test.corr[,i], y = pro1.test[,i])) +
    theme_test() +
    geom_point(colour = "#5B9BD5") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("ADT (Predicted)") +
    ylab("ADT") +
    ggtitle(paste(colnames(rna1)[i], " | ","rho =", round(corrs.sc[i,2], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Single Cell Scatter Plots/", "Corrected.scCF.", i,".png"), dpi = "retina")
}