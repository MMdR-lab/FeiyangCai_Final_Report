library(dplyr)
library(devtools)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(CellChat)
library(jsonlite)

sample <- Load10X_Spatial(data.dir = "D:/McGill/phd course/single-cell/final/FFPE_Human_Skin_Melanoma",
                            filename = "filtered_feature_bc_matrix.h5")
SpatialFeaturePlot(sample, features = "nCount_Spatial") + theme(legend.position = "right")

sample <- SCTransform(sample, assay = "Spatial", verbose = T)

sample <- RunPCA(sample, assay = "SCT", npcs = 50)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:50)
sample <- FindClusters(sample, resolution = 1)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:50)

DimPlot(sample, reduction = "umap", pt.size = 1)

levels(Idents(sample)) <- c("Macrophage", "Tumor PEX3-", "Tumor PEX3-", "Macrophage",
                            "Macrophage", "Macrophage", "T.CD8", "Tumor PEX3+",
                            "Macrophage", "Tumor PEX3-", "T.CD8", "Tumor PEX3-",
                            "Macrophage", "Tumor PEX3+", "Tumor PEX3-", "Tumor PEX3+",
                            "Tumor PEX3+", "Tumor PEX3+", "Macrophage")

FeaturePlot(sample, features = c("MITF", "PEX3", "VEGFB", "THBS3", "TIGIT", "CSF1R",
                                 "CD8A", "CD274", "PVR", "PGF", "CD47", "FLT1"),
            ncol = 3, pt.size = 0.1)

VlnPlot(sample, features = c("MITF", "PEX3", "VEGFB", "THBS3", "TIGIT", "CSF1R",
                             "CD8A", "CD274", "PVR", "PGF", "CD47", "FLT1"),
        ncol = 3, pt.size = 0)

SpatialFeaturePlot(sample, features = c("MITF", "PEX3", "VEGFB", "THBS3", "TIGIT", "CSF1R",
                                        "CD8A", "CD274", "PVR", "PGF", "CD47", "FLT1"),
                   ncol = 4, alpha = c(0.1, 1))
SpatialDimPlot(sample, label = F, label.size = 3, alpha = c(1, 1))
SpatialDimPlot(sample, label = F, label.size = 3, alpha = c(1, 1), group.by = "SCT_snn_res.1")


saveRDS(sample, "D:/McGill/phd course/single-cell/final/FFPE_Human_Skin_Melanoma/sample.rds")

# Cellchat v2
data.input <- GetAssayData(sample, slot = "data", assay = "SCT")
meta <- data.frame(labels = Idents(sample), row.names = names(Idents(sample)))
spatial.locs <- GetTissueCoordinates(sample, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs <- spatial.locs[, -3]
scale.factors <- fromJSON(txt = file.path("D:/McGill/phd course/single-cell/final/FFPE_Human_Skin_Melanoma/spatial",
                                          "scalefactors_json.json"))
scale.factors <- list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres,
                      fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef,
                      lowres = scale.factors$tissue_lowres_scalef)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = T,
                              scale.distance = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "D:/McGill/phd course/single-cell/final/FFPE_Human_Skin_Melanoma/output/cellchat.rds")

groupSize <- as.numeric(table(cellchat@idents))
netVisual_aggregate(cellchat, signaling = "TIGIT", layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_bubble(cellchat, signaling = "PVR", remove.isolate = F)

netVisual_aggregate(cellchat, signaling = "PD-L1", layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

netVisual_aggregate(cellchat, signaling = "THBS", targets.use = 1, layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = "VEGF", targets.use = 1, layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

netVisual_aggregate(cellchat, signaling = "IL10", layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
