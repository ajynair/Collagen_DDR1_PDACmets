# 2022 Jan 05 created
## script analyses the PDAC primary datasets from 2021 Lee Pancreatic cancer data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156405


setwd("D:/pCloud Sync/Columbia/SchwabeLab/Collaborations/KarinLab/DDR1_PDAC/code/Collagen_DDR1_PDACmets")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


markers <- unlist(strsplit("EPCAM, KRT8, CD3D, IL7R, TOP2A, MKI67, NKG7, GNLY, FCGR3A, CD14, LYZ, CD79A, MS4A1, JCHAIN, HBA1, HBA2, TSPAN13, GPR183, PECAM1, CDH5, ACTA2, COL1A1",split = ", "))
markerGenes <- toupper(c("Cpa1","Ctrb1",'Cela3b', #Acinar
                 'Epcam', 'Cldn3', #Epithelial-Ductal
                 "Gcg","Ttr", #alpha
                 "Ins1","Ins2", #beta
                 # "Sst","Lepr", #delta
                 "Ppy", #"Stmn2", # gamma
                 # "Ghrl", # epsilon
                 "Cfd","Car3","Cyp2e1", #Epithelial metaplasia #cluster 25
                 "Pdgfra","Col1a1","Col3a1", #Stellate cells
                 'Msln', 'Upk3b', #Mesothelial
                 "Pecam1","Kdr", # Endothelial
                 "Mmrn1", "Pdpn", #lymphatic Endothelial
                 
                 "Ctss","Lyz2", # Macrphages
                 "Cd14","Lgals3", # Monocytes
                 
                 "Cd3d","Il7r", # T cell
                 "Nkg7", # NK cell
                 "Ms4a1","Cd79a", # B cell
                 "Mzb1","Jchain", # Plasma
                 "Hba-a1","Snca", # Erythroid
                 
                 "Top2a","Mki67" ) )
#******************************************************#******************************************************
#****************************************************** read the primary datasets
#******************************************************#******************************************************
loc <- "data/2021_Lee"
dir <- list.dirs(path = loc,full.names = F, recursive = F)
dir
dir <- dir[grep("^P",dir)] #keep only the primary samples
dir
i=dir[1]
seurat.list <- vector(mode = "list", length = length(dir)) 
names(seurat.list) <- dir
seurat.list
# name.list <- rep(NA, length(dir))
i=1
for (i in 1:length(dir)) {
  if (file.exists(file.path(loc, dir[i]))){
    sampleName <- substr(dir[i],1,2)
    
    # name.list[i] <- sampleName
    print(sampleName)
    dat_filt <- Read10X(file.path(loc, dir[i]))
    dim(dat_filt)
    
    seurat <- CreateSeuratObject(counts = dat_filt, project = sampleName, min.cells = 3, min.features = 200)
    seurat
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RPL|^RPS|MRPS|^MRPL")#
    seurat <- subset(seurat, subset = percent.mt < 30 & nCount_RNA < 30000 )#& nFeature_RNA > 200 & nFeature_RNA < 6500 )#&
    seurat.list[i] <- seurat
    # names(seurat.list[i]) <- sampleName
  }
}
length(seurat.list)
seurat.list[sapply(seurat.list, is.null)] <- NULL
seurat.list


# *************************************************************
# combining all the samples ; mito threshold
# *************************************************************
tmp <- names(seurat.list)
tmp 
# required samples
req_samples <-tmp #[-which(tmp %in% c("S41","S42","S43","S44"))] 
req_sample_pos <- which(tmp %in% req_samples)


req.list <- seurat.list[req_sample_pos]
names(req.list)

req.list <- lapply(X = req.list, FUN = function(x) {
  print(unique(x@meta.data$orig.ident))
  # x <- subset(x, subset = percent.mt < 40 )#& nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA < 40000 )#&
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE, nfeatures = 3000)
})


names(req.list)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = req.list, nfeatures = 3000)

# scale and compute pca for RPCA integration method
req.list <- lapply(X = req.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# perform RPCA integration of all samples
date(); data.anchors <- FindIntegrationAnchors(object.list = req.list, reduction = "rpca", anchor.features = features); date(); # reference = ref_sample_pos
# create an 'integrated' data assay
data.combined_all <- IntegrateData(anchorset = data.anchors) #, dims = 1:50
date()


DefaultAssay(data.combined_all) <- "integrated"
# Run the standard workflow for visualization and clustering
data.combined_all <- ScaleData(data.combined_all, verbose = FALSE)
data.combined_all <- RunPCA(data.combined_all, npcs = 50, verbose = FALSE)
ElbowPlot(data.combined_all, ndims = 50)
ndim <- 10
data.combined_all <- RunUMAP(data.combined_all, reduction = "pca", dims = 1:ndim)
data.combined_all <- FindNeighbors(data.combined_all, reduction = "pca", dims = 1:ndim)
data.combined_all <- FindClusters(data.combined_all)#, resolution = 0.5
DimPlot(data.combined_all, label = T)
DimPlot(data.combined_all,group.by = "orig.ident", label = F,shuffle = T)
FeaturePlot(object = data.combined_all, features = c("percent.mt"))
FeaturePlot(object = data.combined_all, features = c("percent.ribo"))
FeaturePlot(object = data.combined_all, features = c("nCount_RNA"))
FeaturePlot(object = data.combined_all, features = c("nFeature_RNA"))
VlnPlot(object = data.combined_all, features = c("percent.mt"))
VlnPlot(object = data.combined_all, features = c("percent.ribo"))
VlnPlot(object = data.combined_all, features = c("nCount_RNA"))
VlnPlot(object = data.combined_all, features = c("nFeature_RNA"))

DefaultAssay(data.combined_all) <- "RNA"
DotPlot(data.combined_all, features = markers) + RotatedAxis()
DotPlot(data.combined_all, features = markers, cluster.idents = T) + RotatedAxis()


# get cluster markers
table(data.combined_all@active.ident)
markers <- FindAllMarkers(data.combined_all,only.pos = T, logfc.threshold = 0.25)#,max.cells.per.ident = 1000
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
head(markers,20)
# write.csv(markers, file = "results/markers_pdac_primary_clusters.csv")


# giving cell types
DimPlot(data.combined_all,label = T)
# data.combined_all[["Cell.clusters"]] <- Idents(object = data.combined_all)
DimPlot(data.combined_all,label = T,group.by = "Cell.clusters")
data.combined_all <- SetIdent(data.combined_all, value = "Cell.clusters")
data.combined_all <- RenameIdents(object = data.combined_all, 
                                  '1' = "Epithelial-Tumor",
                                  '4' = "Epithelial-Tumor",
                                  '5' = "Epithelial-Tumor",
                                  '6' = "Epithelial-Tumor",
                                  '10' = "Epithelial-Tumor",
                                  '11' = "Epithelial-Tumor",
                                  '15' = "Epithelial-Tumor",
                                  '16' = "Epithelial-Tumor",
                                  '19' = "Epithelial-Tumor",
                                  
                                  '20' = "Epithelial",
                                  
                                  '17' = "Fibroblasts",
                                  
                                  '2' = "Myeloid",
                                  '9' = "Myeloid",
                                  '13' = "Myeloid",
                                  
                                  '14' = "B_Plasma",
                                  
                                  '0' = "T_NK",
                                  '3' = "T_NK",
                                  '8' = "T_NK",
                                  '7' = "T_NK",
                                  
                                  '18' = "Mixed",
                                  '12' = "Mixed")


data.combined_all[["Cell.types"]] <- Idents(object = data.combined_all)
DimPlot(data.combined_all,label = T, group.by = "Cell.types")
table(data.combined_all@meta.data$Cell.types)
table(data.combined_all@meta.data$Cell.clusters)
table(data.combined_all@meta.data$cell.type)

DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters", label = F)

DotPlot(data.combined_all, features = markerGenes) + RotatedAxis()

# saveRDS(data.combined_all, file = "results/seuratObj_pdac_primary.rds")


# **************************************** M1/M2 analysis

m1Sig <- toupper(unlist(strsplit("Azin1, Cd38, Cxcl10, Cxcl9, Fpr2, Il18, Il1b, Irf5, Nifkbiz, Tlr4, Tnf, Cd80",split=", ")))
m2Sig <- toupper(unlist(strsplit("Alox5, Arg1, Chil3, Cd163, Il10, Il10ra, Il10rb, Irf4, Kif4, Mrc1, Myc, Socs2, Tgm2",split=", ")))

cells = WhichCells(object = data.combined_all, idents = "Myeloid")

genes <- m1Sig
genes <- intersect(genes,rownames(data.combined_all@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = data.combined_all, features = c(gene), max.cutoff = "q99",pt.size = 1, order = T)


DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters", label = T,cells = cells)

sig <- m1Sig
sigName <- "m1Sig"
callName <- paste0(sigName,"1")
data.combined_all <- AddModuleScore(data.combined_all,features = list(sig), name = sigName)

FeaturePlot(data.combined_all, callName, max.cutoff = "q95",pt.size = 1,order = T,cells = cells) + scale_color_viridis_c(option = "viridis") +labs(title = "Human PDAC primary M1 signature")
VlnPlot(data.combined_all, group.by = "Cell.clusters", features = callName, idents = "Myeloid")


genes <- m2Sig
genes <- intersect(genes,rownames(data.combined_all@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = data.combined_all, features = c(gene), max.cutoff = "q99",pt.size = 1, order = T)

sig <- m2Sig
sigName <- "m2Sig"
callName <- paste0(sigName,"1")
data.combined_all <- AddModuleScore(data.combined_all,features = list(sig), name = sigName)


FeaturePlot(data.combined_all, callName, max.cutoff = "q99",pt.size = 1,order = T,cells = cells) + scale_color_viridis_c(option = "viridis") +labs(title = "Human PDAC primary M2 signature")
VlnPlot(data.combined_all, group.by = "Cell.clusters", features = callName, idents = "Myeloid")


DimPlot(data.combined_all,label = T,group.by = "Cell.clusters")
data.combined_all <- SetIdent(data.combined_all, value = "Cell.clusters")
data.combined_all <- RenameIdents(object = data.combined_all, 
                                  '1' = "Epithelial-Tumor",
                                  '4' = "Epithelial-Tumor",
                                  '5' = "Epithelial-Tumor",
                                  '6' = "Epithelial-Tumor",
                                  '10' = "Epithelial-Tumor",
                                  '11' = "Epithelial-Tumor",
                                  '15' = "Epithelial-Tumor",
                                  '16' = "Epithelial-Tumor",
                                  '19' = "Epithelial-Tumor",
                                  
                                  '20' = "Epithelial",
                                  
                                  '17' = "Fibroblasts",
                                  
                                  '2' = "M2-like",
                                  '9' = "M1-like",
                                  '13' = "M1-like",
                                  
                                  '14' = "B_Plasma",
                                  
                                  '0' = "T_NK",
                                  '3' = "T_NK",
                                  '8' = "T_NK",
                                  '7' = "T_NK",
                                  
                                  '18' = "Mixed",
                                  '12' = "Mixed")


data.combined_all[["Cell.types.M1M2"]] <- Idents(object = data.combined_all)
DimPlot(data.combined_all,label = T, group.by = "Cell.types.M1M2")


# saveRDS(data.combined_all, file = "results/seuratObj_pdac_primary.rds")

# getting the cell type annotation
colnames(data.combined_all@meta.data)
head(data.combined_all@meta.data)
celltypes <- data.combined_all@meta.data[,c(1,2,3,4,8,9,10)]
head(celltypes)
# write.csv(celltypes, file = "results/cellAnnotation_pdac_primary.csv")


# *************************************************************
# using saved data
# *************************************************************
data.combined_all <- readRDS(file = "results/seuratObj_pdac_primary.rds")
DimPlot(data.combined_all, label = T)
DimPlot(data.combined_all,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types.M1M2", label = F)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters", label = F)

DefaultAssay(data.combined_all) <- "RNA"
DotPlot(data.combined_all, features = markerGenes) + RotatedAxis()
DotPlot(data.combined_all, features = markerGenes, cluster.idents = T) + RotatedAxis()


data.combined_all <- SetIdent(data.combined_all, value = "Cell.types")
req_subset <- subset(data.combined_all, idents = c("Mixed"),invert=T)
req_subset
DimPlot(req_subset,label = T, group.by = "ident")

#******************************************************#******************************************************
#******************************************************PLOTS 
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")


table(data.combined_all@active.ident)
data.combined_all <- SetIdent(data.combined_all, value = "Cell.types")
req_subset <- subset(data.combined_all, idents = c("Mixed"),invert=T)
req_subset
DimPlot(req_subset,label = T, group.by = "ident")#,pt.size = 2
DimPlot(req_subset,group.by = "orig.ident", label = F,shuffle = T,)


# CELL TYPES UMAP
data_req <- req_subset
DimPlot(object = data_req, group.by = "ident", label = T)
levels(data_req@active.ident)


# cols <- c("#00468BFF", "#ADB6B6FF",   "#0099B4FF", "#1B1919FF", "#ED0000FF",      "#AD002AFF",    "#42B540FF",     "#FF7F00")

cols <- c("#FF0000FF","#CC33FFFF","#1B1919FF","#FDBF6F","#42B540FF","#AD002AFF")
DimPlot(object = req_subset, group.by = "Cell.types", label = T, cols = cols)

DimPlot(object = req_subset, group.by = "Cell.types", label = T, cols = cols) +  labs(title = "Human PDAC primary")
ggsave(filename = paste0("results/figs/hum_PDACprimary_celltypes_label.pdf"))

DimPlot(object = req_subset, group.by = "Cell.types", label = F, cols = cols) +   NoLegend()  + labs(title = "Human PDAC primary")
ggsave(filename = paste0("results/figs/hum_PDACprimary_celltypes_nolegend.pdf"))
# ggsave(filename = paste0("results/figs/hum_PDACprimary_celltypes_legend.eps"))

DimPlot(object = req_subset, group.by = "orig.ident", label = F, shuffle = T) +   scale_color_aaas()  + labs(title = "Human PDAC primary")
ggsave(filename = paste0("results/figs/hum_PDACprimary_samples.pdf"))


DefaultAssay(req_subset) <- "RNA"

# dOTPLOTS
markers <- unlist(strsplit("EPCAM, KRT8, CPA1, CTRB1, ACTA2, COL1A1, FCGR3A, CD14, LYZ, CD79A, MS4A1, JCHAIN, CD3D, IL7R, NKG7, GNLY",split = ", "))
DotPlot(req_subset, features = markers) + RotatedAxis() & scale_color_distiller(palette = "RdBu")
ggsave(filename = "results/figs/hum_PDACprimary_dotplot_markers.pdf", width = 7, height = 3.8)




# genes <- unlist(strsplit("EPCAM, KRT19, DDR1, LAIR1, MRC2, ITGB1, CD24, CD44", split = ", "))
genes <- unlist(strsplit("EPCAM, KRT19, DDR1, DDR2, LAIR1, MRC2, ITGB1, NFE2L2, SDC1, TFAM, CDC42", split = ", "))
genes
genes <- intersect(genes,rownames(req_subset@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q99",pt.size = 1, order = F)

for (gene in genes) {
  FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q95",pt.size = 1, order = F) #+ scale_color_viridis_c()
  ggsave(filename = paste0("results/figs/hum_PDACprimary_gene_",gene,".pdf"))
}



genes <- unlist(strsplit("MMP1,MMP2,MMP8,MMP9,MMP13,MMP14", split = ","))
genes <- intersect(genes,rownames(req_subset@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q99",pt.size = 1, order = F)

for (gene in genes) {
  FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q95",pt.size = 1, order = F) #+ scale_color_viridis_c()
  ggsave(filename = paste0("results/figs/hum_PDACprimary_gene_",gene,".pdf"))
}


# **************************************** M1/M2
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.M1M2", label = F)

req_subset <- SetIdent(req_subset, value = "Cell.types") 
cells = WhichCells(object = req_subset, idents = "Myeloid")


DimPlot(object = req_subset, group.by = "Cell.types", label = F, cols = cols[4], cells = cells, pt.size = 3)+ xlim(-2,8) +ylim(5,15)
ggsave(filename = paste0("results/figs/hum_PDACprimary_M1M2_celltypes.pdf"))

req_subset <- SetIdent(req_subset, value = "Cell.types.M1M2")
col1 <- c("#FF9900FF", "#6699FFFF")
DimPlot(object = req_subset, group.by = "Cell.types.M1M2", label = F, cols = col1, cells = cells, pt.size = 3)+ xlim(-2,8) +ylim(5,15)
ggsave(filename = paste0("results/figs/hum_PDACprimary_M1M2_celltypes_M1M2.pdf"))


req_subset <- SetIdent(req_subset, value = "Cell.types") 
sigName <- "m1Sig"
callName <- paste0(sigName,"1")
VlnPlot(req_subset, group.by = "Cell.types.M1M2", features = callName, idents = "Myeloid",cols = col1)
ggsave(filename = paste0("results/figs/hum_PDACprimary_M1M2_M1sig.pdf"), width = 5, height = 8)

sigName <- "m2Sig"
callName <- paste0(sigName,"1")
VlnPlot(req_subset, group.by = "Cell.types.M1M2", features = callName, idents = "Myeloid",cols = col1)
ggsave(filename = paste0("results/figs/hum_PDACprimary_M1M2_M2sig.pdf"), width = 5, height = 8)


genes <- unlist(strsplit("MMP1,MMP2,MMP8,MMP9,MMP13,MMP14", split = ","))
genes <- intersect(genes,rownames(req_subset@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q99", order = T, cells = cells, pt.size = 3)+ xlim(-2,8) +ylim(5,15)

for (gene in genes) {
  FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q95", order = T, cells = cells, pt.size = 3)+ xlim(-2,8) +ylim(5,15)
  ggsave(filename = paste0("results/figs/hum_PDACprimary_M1M2_gene_",gene,".pdf"))
}

