BiocManager::install("monocle")

library("Seurat")
library("monocle")
library("CellChat")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("msigdbr")
library("fgsea")
library("GSVA")
library("GSEABase")

melanoma.data <- read.csv("D:/McGill/phd course/single-cell/final/GSE115978/GSE115978_counts.csv")
rownames(melanoma.data) <- melanoma.data$X
melanoma.data <- melanoma.data[, -1]
metadata <- read.csv("D:/McGill/projects/project BRAFi/lighthouse/GSE115978/GSE115978_cell.annotations.csv")

melanoma <- CreateSeuratObject(counts = melanoma.data, project = "melanoma")
melanoma$samples <- metadata$samples
melanoma$cell.type <- metadata$cell.types
melanoma$treatment <- metadata$treatment.group

melanoma <- SCTransform(melanoma, assay = "RNA", verbose = T)
saveRDS(melanoma, "D:/McGill/projects/project BRAFi/lighthouse/GSE115978/sample/melanoma.rds")

melanoma <- readRDS("D:/McGill/projects/project BRAFi/lighthouse/GSE115978/sample/melanoma.rds")

unique(melanoma$cell.type)

melanoma <- NormalizeData(melanoma, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
melanoma <- FindVariableFeatures(melanoma, selection.method = "vst", # nfeatures = 200,
                                 mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(melanoma), 10)

#Scaling the data standardization
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(melanoma)
melanoma <- ScaleData(melanoma, features = all.genes)

melanoma <- RunPCA(melanoma, npcs = 50, features = VariableFeatures(object = melanoma))
# Examine and visualize PCA results a few different ways
print(melanoma[["pca"]], dims = 1:5, nfeatures = 5)

melanoma <- FindNeighbors(melanoma, dims = 1 : 50)
melanoma <- FindClusters(melanoma, resolution = 0.5)


# UMAP
melanoma = RunUMAP(object = melanoma, reduction = "pca", dims = 1 : 50)
DimPlot(melanoma, reduction = "umap", group.by = "cell.type")
DimPlot(melanoma, reduction = "umap", group.by = "samples")

TME <- subset(melanoma, subset = cell.type %in%
                c("B.cell", "CAF", "Endo.", "Macrophage", "NK", "T.CD4", "T.CD8", "T.cell"))
tumor <- subset(melanoma, subset = cell.type == "Mal")

DimPlot(TME, reduction = "umap", group.by = "cell.type")
DimPlot(TME, reduction = "umap", group.by = "samples")
DimPlot(tumor, reduction = "umap", group.by = "samples")

remotes::install_github("immunogenomics/harmony")
library(harmony)

# TME
TME <- RunHarmony(TME, group.by.vars = "samples", project.dim = FALSE)
TME <- FindNeighbors(TME, dims = 1 : 50, reduction = "harmony")
TME <- FindClusters(TME, resolution = 0.5)

# UMAP
TME = RunUMAP(object = TME, dims = 1 : 50, reduction = "harmony")
TME = subset(TME, idents = c(0:10, 13))
DimPlot(TME, reduction = "umap", group.by = "cell.type")
DimPlot(TME, reduction = "umap")
DimPlot(TME, reduction = "umap", group.by = "samples")

FeaturePlot(TME, features = c("CD19", "PDGFRA", "PECAM1", "CSF1R",
                              "CCR7", "CD4", "CD8A", "GZMB", "FOXP3", "MKI67",
                              "FCGR3A"),
            ncol = 3, pt.size = 0.1)

markers <- FindMarkers(TME, ident.1 = 13, ident.2 = 2)

levels(Idents(TME)) <- c("T.CD8", "T.naive", "B.cell", "Macrophage", "T.CD4", "T.CD4.reg",
                         "NK", "T.CD8.mitotic", "CAF", "T.CD8", "Endothelial", "B.cell")

saveRDS(TME, "D:/McGill/phd course/single-cell/final/GSE115978/sample/TME.rds")

TME_pre <- subset(TME, subset = treatment == "treatment.naive")
TME_post <- subset(TME, subset = treatment == "post.treatment")

DimPlot(TME_pre, reduction = "umap")
DimPlot(TME_post, reduction = "umap")

# CD8 T cell
CD8 <- subset(TME, idents = "T.CD8")

CD8 <- FindNeighbors(CD8, dims = 1 : 30, reduction = "harmony")
CD8 <- FindClusters(CD8, resolution = 0.5)

# UMAP
CD8 = RunUMAP(object = CD8, dims = 1 : 30, reduction = "harmony")
DimPlot(CD8, reduction = "umap")

Exhaustion_score <- apply(CD8@assays[["SCT"]]@data[c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "BATF", "IRF4", "TNFRSF9"),], 2, mean)
CD8$Exhaustion_score <- Exhaustion_score

CD8_pre <- subset(CD8, subset = treatment == "treatment.naive")
CD8_post <- subset(CD8, subset = treatment == "post.treatment")

DimPlot(CD8_pre, reduction = "umap")
DimPlot(CD8_post, reduction = "umap")

FeaturePlot(CD8_pre, features = "Exhaustion_score",
            ncol = 2, pt.size = 1, max.cutoff = 2.2)
FeaturePlot(CD8_post, features = "Exhaustion_score",
            ncol = 2, pt.size = 1, max.cutoff = 2.2)


FeaturePlot(CD8_pre, features = c("PDCD1", "HAVCR2"),
            ncol = 2, pt.size = 1, max.cutoff = 3.9)
FeaturePlot(CD8_post, features = c("PDCD1", "HAVCR2"),
            ncol = 2, pt.size = 1, max.cutoff = 3.9)

FeaturePlot(CD8_pre, features = c("LAG3", "BATF"),
            ncol = 2, pt.size = 1, max.cutoff = 2.5)
FeaturePlot(CD8_post, features = c("LAG3", "BATF"),
            ncol = 2, pt.size = 1, max.cutoff = 2.5)

FeaturePlot(CD8_pre, features = c("TIGIT"), pt.size = 1, max.cutoff = 4.3)
FeaturePlot(CD8_post, features = c("TIGIT"), pt.size = 1, max.cutoff = 4.3)

FeaturePlot(CD8_pre, features = c("GZMB", "PRF1"),
            ncol = 2, pt.size = 1, max.cutoff = 4.5)
FeaturePlot(CD8_post, features = c("GZMB", "PRF1"),
            ncol = 2, pt.size = 1, max.cutoff = 4.5)

FeaturePlot(CD8_pre, features = c("IRF4", "TNFRSF9"),
            ncol = 2, pt.size = 1, max.cutoff = 3.9)
FeaturePlot(CD8_post, features = c("IRF4", "TNFRSF9"),
            ncol = 2, pt.size = 1, max.cutoff = 3.9)

# macrophage
macrophage <- subset(TME, idents = "Macrophage")

macrophage <- FindNeighbors(macrophage, dims = 1 : 30, reduction = "harmony")
macrophage <- FindClusters(macrophage, resolution = 0.1)

# UMAP
macrophage = RunUMAP(object = macrophage, dims = 1 : 30, reduction = "harmony")
DimPlot(macrophage, reduction = "umap")

FeaturePlot(macrophage, features = "CD274", pt.size = 1)
macrophage$treatment <- factor(macrophage$treatment, levels = c("treatment.naive", "post.treatment"))
VlnPlot(macrophage, features = "CD274", group.by = "treatment", pt.size = 0)

CD8_pre <- subset(CD8, subset = treatment == "treatment.naive")
CD8_post <- subset(CD8, subset = treatment == "post.treatment")






# tumor
tumor <- RunHarmony(tumor, group.by.vars = "samples", project.dim = FALSE, dims = 1:50)
tumor <- FindNeighbors(tumor, dims = 1 : 50, reduction = "harmony")
tumor <- FindClusters(tumor, resolution = 0.2)


# UMAP
tumor = RunUMAP(object = tumor, reduction = "harmony", dims = 1 : 50)
tumor <- subset(tumor, idents = c(0:4, 6))
DimPlot(tumor, reduction = "umap", pt.size = 1)
DimPlot(tumor, reduction = "umap", group.by = "samples", pt.size = 1)
DimPlot(tumor, reduction = "umap", group.by = "treatment", pt.size = 1)


markers <- FindMarkers(tumor, ident.1 = "RXRG pigmented", ident.2 = "pigmented")
write.csv(markers, "D:/McGill/phd course/single-cell/final/GSE115978/sample/tumor_markers.csv")

FeaturePlot(tumor, features = c("PMEL", "MLANA", "CCL5", "FGF5",
                                "CDK1", "MKI67", "JUN", "FOS",
                                "FOSB", "MITF", "RXRG"),
            ncol = 3)
levels(Idents(tumor)) <- c("pigmented", "immune", "proliferative", "pigmented", "AP-1", "RXRG pigmented")

FeaturePlot(tumor, features = "MIF", pt.size = 1)
VlnPlot(tumor, features = "MIF", group.by = "pseudo_state", pt.size = 0)

saveRDS(tumor, "D:/McGill/phd course/single-cell/final/GSE115978/sample/tumor.rds")

# Pseudotime
Mono_matrix <- as(as.matrix(GetAssayData(tumor, slot = "counts")), 'sparseMatrix')
feature_ann <- data.frame(gene_id = rownames(Mono_matrix),
                          gene_short_name = rownames(Mono_matrix))
rownames(feature_ann) <- rownames(Mono_matrix)
Mono_fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- tumor@meta.data
rownames(sample_ann) <- colnames(Mono_matrix)
Mono_pd <- new("AnnotatedDataFrame", data = sample_ann)
Mono.cds <- newCellDataSet(Mono_matrix, phenoData = Mono_pd, featureData = Mono_fd,
                           expressionFamily = negbinomial.size())
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)

signature <- read.csv("D:/McGill/phd course/single-cell/final/signature.csv")
ordering_genes <- c(signature$pigmentation, signature$invasion, signature$neuro, signature$MITF.targets)
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Mono.cds <- setOrderingFilter(Mono.cds, intersect(unsup_clustering_genes$gene_id, ordering_genes))
Mono.cds <- reduceDimension(Mono.cds, max_components = 2, method = 'DDRTree')
Mono.cds <- orderCells(Mono.cds)

plot_cell_trajectory(Mono.cds, cell_size = 1)
plot_cell_trajectory(Mono.cds, cell_size = 1) +
  facet_wrap("~treatment", nrow = 1)

plot_cell_trajectory(Mono.cds, color_by = "Pseudotime", size = 1, show_backbone = F)
plot_cell_trajectory(Mono.cds, color_by = "Pseudotime", size = 1, show_backbone = T) +
  facet_wrap("~treatment", nrow = 1)


Mono.cds$PEX3 <- tumor@assays[["SCT"]]@data["PEX3",]
plot_cell_trajectory(Mono.cds, color_by = "PEX3", cell_size = 0.5, show_backbone = T) +
  scale_color_gradient(low = "grey", high = "blue")

levels(Mono.cds$State) <- c(rep("Pigmented", 6), "Stem-like")
PEX3_state <- as.character(Mono.cds$State)
PEX3_pos <- intersect(which(tumor@assays[["SCT"]]@data["PEX3",] > 0),
                      which(PEX3_state == "Pigmented"))
PEX3_state[PEX3_pos] <- rep("Pigmented_PEX3+", length(PEX3_pos))
PEX3_pos <- intersect(which(tumor@assays[["SCT"]]@data["PEX3",] > 0),
                      which(PEX3_state == "Stem-like"))
PEX3_state[PEX3_pos] <- rep("Stem-like_PEX+", length(PEX3_pos))
Mono.cds$State <- as.factor(PEX3_state)

plot_genes_in_pseudotime(Mono.cds[c("PMEL", "NGFR", "PEX3"),], color_by = "State")
plot_genes_in_pseudotime(Mono.cds[c("AGPS", "PEX19", "PEX16", "PEX3",
                                    "ABCD3", "SCP2", "PEX1", "PEX2",
                                    "PEX13", "FAR1", "ATG12", "DNM1L"),], color_by = "State")


saveRDS(Mono.cds, "D:/McGill/phd course/single-cell/final/GSE115978/output/pseudotime_tumor.rds")

plot_pseudotime_heatmap(Mono.cds[c("PEX3", intersect(unsup_clustering_genes$gene_id, ordering_genes)),],
                        num_clusters = 4, show_rownames = T, return_heatmap = T)

tumor$pseudo_state <- Mono.cds$State

# CellChat
melanoma.data <- merge(tumor, TME)
exprs_mtx <- melanoma.data@assays[["SCT"]]@data
metadata <- as.data.frame(c(tumor$pseudo_state, Idents(TME)))

cellchat <- createCellChat(object = exprs_mtx)
cellchat <- addMeta(cellchat, meta = metadata, meta.name = "cellType")
cellchat <- setIdent(cellchat, ident.use = "cellType")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(5, 6, 8:12), remove.isolate = F)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5, 6, 8:12), remove.isolate = F)

netVisual_aggregate(cellchat, signaling = "MHC-I", vertex.weight = groupSize)
netVisual_bubble(cellchat, sources.use = c(1, 2), signaling = "MHC-I", remove.isolate = F)

netVisual_aggregate(cellchat, sources.use = c(1:4), signaling = "CCL", vertex.weight = groupSize)
netVisual_bubble(cellchat, sources.use = c(1:4), signaling = "CCL", targets.use = 8, remove.isolate = F)


netVisual_aggregate(cellchat, sources.use = c(1:4), signaling = "CXCL", vertex.weight = groupSize)


netVisual_aggregate(cellchat, signaling = "MIF", vertex.weight = groupSize)
netVisual_bubble(cellchat, sources.use = c(1:4), signaling = "MIF", targets.use = 8, remove.isolate = F)

netVisual_aggregate(cellchat, signaling = "PD-L1", vertex.weight = groupSize)
netVisual_bubble(cellchat, sources.use = 8, signaling = "PD-L1", remove.isolate = F)

netVisual_aggregate(cellchat, signaling = "TIGIT", vertex.weight = groupSize)
netVisual_bubble(cellchat, sources.use = c(5, 12), signaling = "TIGIT", targets.use = c(1:4), remove.isolate = F)

netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = 8, remove.isolate = F)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = 8,  signaling = "THBS", remove.isolate = F)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = 8,  signaling = "IL10", remove.isolate = F)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = 8,  signaling = "VEGF", remove.isolate = F)

netVisual_bubble(cellchat, sources.use = c(5,12), targets.use = c(1:4), remove.isolate = F)

saveRDS(cellchat, "D:/McGill/phd course/single-cell/final/GSE115978/output/cellchat.rds")

# tumor GSEA
Idents(tumor) <- tumor$pseudo_state
markers <- FindAllMarkers(tumor, only.pos = F, logfc.threshold = 1)
write.csv(markers, "D:/McGill/phd course/single-cell/final/GSE115978/output/GSEA/markers.csv")

d <- markers[which(markers$cluster == "Pigmented"),]
d <- markers[which(markers$cluster == "Pigmented_PEX3+"),]
d <- markers[which(markers$cluster == "Stem-like"),]
d <- markers[which(markers$cluster == "Stem-like_PEX+"),]

names(d)[1] <- "SYMBOL"
rownames(d) <- d$SYMBOL

id <- bitr(rownames(d), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

geneList <- merge(d, id, by = "SYMBOL", all = FALSE)
geneList = geneList[order(geneList$avg_log2FC, decreasing = TRUE),]

gene.expr <- geneList$avg_log2FC
names(gene.expr) <- geneList$ENTREZID

go <- gseGO(gene.expr, ont = "ALL", OrgDb = org.Hs.eg.db)
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
gseaplot2(go, 
          go@result$ID[28],
          title = go@result$Description[28],
          base_size = 15,
          color = "green",
          pvalue_table = T,
          ES_geom = "line")
write.csv(sortgo, "D:/McGill/phd course/single-cell/final/GSE115978/output/GSEA/Stem-like_PEX3.csv")
write.csv(go@result, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_GSEA.csv")

gene <- go@result[["core_enrichment"]][28]
gene <- unlist(strsplit(gene, split = "/"))
gene <- bitr(gene, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
