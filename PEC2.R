##### Packages #####
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(umap)

##### Functions #####

# Head2: shows a subset of the data.frame/matrix, without showing all the columns
head2 <- function(m, n = 6) m[1:n, 1:n]



# ------------------------------------------------------------------------------ #


##### Obtain Data #####
# Data from publication: DOI 10.15252/msb.20188503

pro <- read.table("Data/Proteomic_Data.csv", sep = ";", header = T, dec = ".")
rna <- read.table("Data/Transcriptomic_Data.csv", sep = ";", header = T)

## RNA and Protein levels distribution

pro.plot <- select(pro, -Gene.name, -Gene.ID) %>% gather
rna.plot <- select(rna, -Gene.name, -Gene.ID) %>% gather

rna.plot$Source <- "RNA"
pro.plot$Source <- "Protein"
data.plot <- rbind(rna.plot, pro.plot)


p1 <- ggplot(pro.plot, aes(x = key, y = value)) +
  theme_bw() +
  geom_boxplot(fill = "#A9D18E") +
  ggtitle("Protein") +
  ylab("log10(iBAQ)") +
  scale_y_continuous(trans = "log10") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))


p2 <- ggplot(rna.plot, aes(x = key, y = value)) +
  theme_bw() +
  geom_boxplot(fill = "#9DC3E6") +
  ggtitle("RNA") +
  ylab("log10(FPKM)") +
  scale_y_continuous(trans = "log10") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))


p3 <- ggplot(data.plot, aes(group = Source, x = value, fill = Source)) +
  theme_bw() +
  geom_density() +
  xlab("log10(Abundance)") +
  ggtitle("") +
  scale_x_continuous(trans = "log10") +
  theme(legend.position = c(.85,.94),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(hjust = .5),
        axis.title.y = element_blank()) +
  scale_color_manual(values = c("Protein" = "#A9D18E", "RNA" = "#9DC3E6"), aesthetics = "fill")

p1.2 <- ggarrange(p1, p2, nrow = 2, labels = LETTERS[1:2])

ggarrange(p1.2, p3, ncol = 2, labels = c("", "C"), widths = c(2,1))

ggsave("RNA & PRO distributions.png",device = "png", dpi = "retina")


##### Integrate Data Sets #####
## Merge Datasets

# Filter genes that are common to both data sets.
# I also removed the column Gene.name from the data.
rna.data <- rna %>% filter(Gene.ID %in% pro$Gene.ID) %>% dplyr::select(-Gene.name)
pro.data <- pro %>% filter(Gene.ID %in% rna$Gene.ID) %>% dplyr::select(-Gene.name)

# Change column names to avoid duplicates
colnames(rna.data)[2:30] <- paste0(colnames(rna.data)[2:30], ".RNA")
colnames(pro.data)[2:30] <- paste0(colnames(pro.data)[2:30], ".PRO")

# Merge data sets
rna.pro <- full_join(rna.data, pro.data)

# Change Gene.ID column to row names
rownames(rna.pro) <- rna.pro$Gene.ID

# Remove Gene.ID (Now in row names) and keep data values
rna.pro <- select(rna.pro, -Gene.ID)


## Remove genes that are present in less than 5 samples 
rna.pro.filt <- apply(rna.pro, 1, function(a) {
  rna <- a[1:29] > 0
  pro <- a[30:58] > 0
  
  # RNA and Protein detected in 5 samples or more
  ifelse(sum(rna & pro) >= 5, return(a), return(NA))

  }) %>% 
  do.call(rbind, .) %>%
  data.frame %>% 
  .[complete.cases(.),]


## Transform values to log scale
rna.pro.filt <- log10(rna.pro.filt + 1)



## PCA

## Vector of gene's position ordered by decreasing variance
rna.max.var <- apply(rna.pro.filt[,01:29], 1, var) %>% order(decreasing = T)
pro.max.var <- apply(rna.pro.filt[,30:58], 1, var) %>% order(decreasing = T)

## PCA of RNA and Protein data
rna.pca <- prcomp(t(rna.pro.filt[rna.max.var[1:1000],01:29]))
pro.pca <- prcomp(t(rna.pro.filt[pro.max.var[1:1000],30:58]))

## Components variance
rna.pca.vars <- (summary(rna.pca)[["sdev"]]^2/sum(summary(rna.pca)[["sdev"]]^2))
pro.pca.vars <- (summary(pro.pca)[["sdev"]]^2/sum(summary(pro.pca)[["sdev"]]^2))

## Create Plots
p1 <- ggplot(data.frame(rna.pca$x), 
       aes(x = PC1, y = PC2,
           label = colnames(rna)[3:31])) +
  theme_bw() +
  geom_point(size = 3, colour = "#9DC3E6") +
  geom_text_repel(max.overlaps = 20) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  xlab(paste0("PC1 (", round(rna.pca.vars[1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(rna.pca.vars[2]*100, 1), "%)")) +
  ggtitle("RNA")

 
p2 <- ggplot(data.frame(pro.pca$x), 
       aes(x = PC1, y = PC2,
           label = colnames(pro)[3:31])) +
  theme_bw() +
  geom_point(size = 3, colour = "#A9D18E") +
  geom_text_repel(max.overlaps = 20) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = .5)) +
  xlab(paste0("PC1 (", round(pro.pca.vars[1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(pro.pca.vars[2]*100, 1), "%)")) +
  ggtitle("Protein")

ggarrange(p1, p2)
ggsave("RNA & PRO PCAs.png",device = "png", dpi = "retina")

## UMAP
xval <- paste0("PC", 1:29) %>% factor(., levels = ., labels = ., ordered = T)

p1.1 <- ggplot(NULL, aes(x = xval, y = rna.pca.vars * 100)) +
  theme_bw() +
  geom_segment(xend = xval, yend = 0) +
  geom_point(size = 3, colour = "#9DC3E6") +
  ylab("% Of Total Variance ") +
  ggtitle("RNA") +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = .5)
  )

rna.UMAP <- umap(rna.pca$x[,1:6])
p1.2 <- ggplot(data.frame(rna.UMAP$layout), 
       aes(x = X1, y = X2,
           label = colnames(rna)[3:31])) +
  theme_bw() +
  geom_point(size = 3, colour = "#9DC3E6") +
  geom_text_repel(max.overlaps = 20) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("RNA")



p2.1 <- ggplot(NULL, aes(x = xval, y = pro.pca.vars * 100)) +
  theme_bw() +
  geom_segment(xend = xval, yend = 0) +
  geom_point(size = 3, colour = "#A9D18E") +
  ylab("% Of Total Variance ") +
  ggtitle("Protein") +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = .5)
  )

pro.UMAP <- umap(pro.pca$x[,1:6])
p2.2 <- ggplot(data.frame(pro.UMAP$layout), 
       aes(x = X1, y = X2,
           label = colnames(rna)[3:31])) +
  theme_bw() +
  geom_point(size = 3, colour = "#A9D18E") +
  geom_text_repel(max.overlaps = 20) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("Protein")

(p1.1 + p2.1)
ggsave("RNA & PRO ScreePlot.png",device = "png", dpi = "retina")

(p1.2 + p2.2)
ggsave("RNA & PRO UMAPS.png",device = "png", dpi = "retina")





##### Protein Prediction from RNA levels #####

## Correction factor
cf <- apply(rna.pro.filt, 1, function(a){
  rna <- a[1:29]
  pro <- a[30:58]
  
  # Flag Cf = Inf and Cf = 0 as NAs
  NAs <- rna == 0 | pro == 0
  
  cf <- pro/rna
  cf[NAs] <- NA
  return(cf)
  }) %>% t

## Protein levels prediction using leave-one-out cross validation method
# Vectors corrs1 and corrs2 keep track of the correlation before and after
# Correction respectively. "rna.corr" stores corrected RNA values for every tissue.

corrs1 <- c()
corrs2 <- c()
rna.corr <- list()

for(i in 1:29){
  # Get RNA and Protein values from tissue i
  t1.rna <- rna.pro.filt[,i]
  t1.pro <- rna.pro.filt[,i+29]

  # Remove 0 values to avoid correlation artifacts
  keep. <- (t1.rna > 0 ) & (t1.pro > 0)
  corrs1[i] <- cor(t1.rna[keep.], t1.pro[keep.])
  
  # Remove tissue i from CFs list
  t1.cf <- cf[,-i]
  
  # Compute median
  t1.med <- apply(t1.cf, 1, function(a) median(a, na.rm = T))
  
  # Predict Protein values from RNA levels using the
  # information from the other 28 tissues.
  t1.rna.corr <- t1.rna * t1.med
  
  # Store values
  rna.corr[[i]] <- t1.rna.corr
  
  # Remove 0 values to avoid correlation artifacts
  keep. <- (t1.rna.corr > 0 ) & (t1.pro > 0)
  corrs2[i] <- cor(t1.rna.corr[keep.], t1.pro[keep.])
}

## Create data set with Protein predicted data
rna.pro.filt.corr <- rna.pro.filt
rna.pro.filt.corr[,1:29] <- do.call(cbind, rna.corr)

## Store correlations before and after correction.
names(corrs1) <- names(corrs2) <- colnames(rna)[3:31]
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

ggsave(filename = paste0("Tissue Scatterplots/", "Boxplot",".png"), dpi = "retina")


# Selection Cell Surface Proteins
# Gene list obtained from: https://doi.org/10.1073/pnas.1808790115
surf.pro <- read.csv("Data/Cell Surface Proteins.csv", header = T, sep = ";")
s.genes <- surf.pro$Ensembl.gene %>% strsplit(., split = ";") %>% purrr::flatten() %>% unlist %>% unique


## Scatter plot before and after correction
tissues <- names(rna)[3:31]

rna.pro.plot <- rna.pro.filt
rna.pro.plot[rna.pro.plot == 0] <- NA

rna.pro.cor.plot <- rna.pro.filt.corr
rna.pro.cor.plot[rna.pro.cor.plot == 0] <- NA

surf.gen <- rownames(rna.pro.filt) %in% s.genes

# Before correction
for (i in 1:29) {
  ggplot(NULL, aes(x = rna.pro.plot[,i], y = rna.pro.plot[,i+29])) +
    theme_test() +
    geom_point(colour = "#7F7F7F") +
    geom_point(aes(x = rna.pro.plot[surf.gen,i], y = rna.pro.plot[surf.gen,i+29]),
               colour = "#2F2F2F") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("RNA") +
    ylab("Protein") +
    ggtitle(paste(tissues[i], " | ","rho =", round(corrs[i,1], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Tissue Scatterplots/", "Uncorrected.", tissues[i],".png"), dpi = "retina")
}

# After Correction
for (i in 1:29) {
  ggplot(NULL, aes(x = rna.pro.cor.plot[,i], y = rna.pro.cor.plot[,i+29])) +
    theme_test() +
    geom_point(colour = "#5B9BD5") +
    geom_point(aes(x = rna.pro.cor.plot[surf.gen,i], y = rna.pro.cor.plot[surf.gen,i+29]),
               colour = "#004F8A") +
    geom_smooth(formula = y~x, method = "lm", 
                se = F, colour = "black", 
                linetype = "dashed", size = .7) +
    xlab("Protein (Predicted)") +
    ylab("Protein") +
    ggtitle(paste(tissues[i], " | ","rho =", round(corrs[i,2], 2))) +
    theme(plot.title = element_text(hjust = .5))
  
  ggsave(filename = paste0("Tissue Scatterplots/", "Corrected.", tissues[i],".png"), dpi = "retina")
}




## Correction factor and Surface proteins
# Compute median Cf for every gene
med <- apply(cf, 1, function(a) median(a, na.rm = T))

# Define groups of increasing Cf value
lim <- quantile(med, seq(0, 1, length.out = 5))
grps <- split(med, cut(med, lim))

# Compute the frequency of surface proteins on every group
surf.genes <- s.genes[s.genes %in% rownames(cf)]
grps.frq <- sapply(grps, function(a) mean(surf.genes %in% names(a)))

# Barplot
grps.plot <- data.frame(Frq = grps.frq, Grps = names(grps.frq))

ggplot(grps.plot, aes(x = Grps, y = Frq*100)) +
  theme_bw() +
  geom_bar(stat = "identity", fill = "#5B9BD5") +
  xlab("Cf") +
  ylab("Surface Proteins (%)") +
  theme(axis.text.x = element_text(size = 10))

ggsave(filename = paste0("Surf Pro barplot",".png"), dpi = "retina")


# Statistical test
values <- sapply(grps, function(a) table(surf.genes %in% names(a)))
chisq.test(values)

## Constant Cf
genes <- c("ENSG00000157500", 
           "ENSG00000110321", 
           "ENSG00000077549", 
           "ENSG00000095139", 
           "ENSG00000196531", 
           "ENSG00000197111")

colnames(cf) <- colnames(rna)[3:31]

cfs.plot <- cf[genes,] %>% data.frame %>% gather("Tissue", "Cf")

gene.names <- rna$Gene.name
names(gene.names) <- rna$Gene.ID
cfs.plot$Gene.ID <- gene.names[genes]

labls <- data.frame(x = 28.5, 
                    y = c(2.67, 2.35, 2.19, 1.83, 1.6, 2),
                    Gene.ID = c("APPL1","ARCN1","CAPZB","EIF4G2","NACA","PCBP2"))


meds <- tapply(cfs.plot$Cf, cfs.plot$Gene.ID, mean)
meds <- data.frame(Gene.ID = names(meds), y = meds)

ggplot(cfs.plot, aes(x = Tissue, y = Cf, group = Gene.ID, colour = Gene.ID)) +
  theme_bw() +
  geom_path(size = 1) +
  geom_hline(data = meds, aes(yintercept = y, colour = Gene.ID),
             linetype = "dashed") +
  geom_text(data = labls, aes(x = x, y = y, label = Gene.ID)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        axis.title.x = element_blank()) +
  ylab("Protein/RNA")

ggsave(filename = paste0("Cf pathplot",".png"), dpi = "retina")