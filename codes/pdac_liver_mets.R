# 2022 Jan 05 created
# script analyses the PDAC Liver-mets dataset from 2021 Lee Pancreatic cancer data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156405


setwd("D:DDR1_PDAC/code/Collagen_DDR1_PDACmets/")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


markerGenes <- c("MYH11","ACTG2","ACTA2", #VSMC
                 "RGS5","RERGL","COLEC11","LRAT", #Hepatic stellate cells
                 "LUM","COL1A1","COL3A1","LOX","TIMP1","DCN", #Myofibroblasts
                 "KDR","AQP1","VWF", #Endothelial
                 "ENG","PECAM1","RAMP3","INMT",#Portal Endo
                 "STAB1",#CV LSEC
                 "SPARCL1","CLEC14A",#PERIPortal  LSEC
                 "TOP2A","MKI67","CENPF", #Cycling
                 'UPK3B','MSLN',  # Mesothelial
                 # 'EPCAM', 'KRT7',"ALB","SERPINA1A",
                 "ALB","APOE","APOB","CPS1", # Hepatocytes
                 "ANXA4","EPCAM", "CFTR","KRT7","KRT8",'KRT19','SOX9', #HCholangiocytes
                 "IL7R","CD8A","CD3D","TRAC","FOXP3", #T cells
                 "NKG7","GNLY", #NK cells
                 "MS4A1","CD79A","IGKC", # B cells
                 "JCHAIN", #Plasma cells
                 "HBB","HBA1","SNCA", #Erythroid/ RBC
                 "CLEC4F","VSIG4","LGALS3", #Kupffer cells
                 "ADGRE1","CD68","CD163","ITGB2","FCGR2A","ITGAM", #Macrophages
                 "FCGR3A","MS4A7", # FCGR3A + Monocytes
                 "CD14","LYZ", # LyZ monocytes
                 "FCER1A", "CST3") # Dendritic cells

#******************************************************#******************************************************
#****************************************************** get the liver met
#******************************************************#******************************************************
data <- Read10X("data/2021_Lee/LiM_filtered_feature_bc_matrix")
dim(data)
sampleName <- "livMets"
# quickly check the data, quality, and cell types
livMets <- CreateSeuratObject(counts = data, project = sampleName, min.cells = 3, min.features = 200)
livMets
pdf(file = paste0("results/reports/sample_",sampleName,".pdf"))
plot.new()
text(x=.5, y=0.5,paste0("Sample: ",data$prefix))
text(x=.5, y=0.4,paste0("Dimensions: ", paste(dim(livMets), collapse = "x")))
livMets[["percent.mt"]] <- PercentageFeatureSet(livMets, pattern = "^MT-")
livMets[["percent.ribo"]] <- PercentageFeatureSet(livMets, pattern = "^RPL|^RPS|MRPS|^MRPL")#
VlnPlot(livMets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(livMets, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(livMets, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(livMets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
livMets <- NormalizeData(object = livMets)
livMets <- FindVariableFeatures(object = livMets, nfeatures = 3000)
livMets <- ScaleData(object = livMets)
livMets <- RunPCA(object = livMets)
ElbowPlot(livMets, ndims = 50)
livMets <- FindNeighbors(object = livMets)
livMets <- FindClusters(object = livMets)
livMets <- RunUMAP(object = livMets, reduction = "pca", dims = 1:30)
FeaturePlot(object = livMets, features = c("percent.mt"))
FeaturePlot(object = livMets, features = c("percent.ribo"))
FeaturePlot(object = livMets, features = c("nCount_RNA"))
FeaturePlot(object = livMets, features = c("nFeature_RNA"))
DimPlot(object = livMets, label = T)
DotPlot(livMets, features = markers) + RotatedAxis()
dev.off()
dev.off()

# make seurat object with proper QC and parameters
livMets <- CreateSeuratObject(counts = data, project = sampleName, min.cells = 3, min.features = 200)
livMets
livMets[["percent.mt"]] <- PercentageFeatureSet(livMets, pattern = "^MT-")
livMets[["percent.ribo"]] <- PercentageFeatureSet(livMets, pattern = "^RPL|^RPS|MRPS|^MRPL")#
livMets <- subset(livMets, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 50000 )#&
VlnPlot(livMets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(livMets, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(livMets, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(livMets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
livMets <- NormalizeData(object = livMets)
livMets <- FindVariableFeatures(object = livMets, nfeatures = 3000)
livMets <- ScaleData(object = livMets)
livMets <- RunPCA(object = livMets)
ElbowPlot(livMets, ndims = 50)
ndims <- 21
livMets <- FindNeighbors(object = livMets, dims = 1:ndims)
livMets <- RunUMAP(object = livMets, reduction = "pca", dims = 1:ndims)
livMets <- FindClusters(object = livMets)

DimPlot(object = livMets, label = T)
DotPlot(livMets, features = markerGenes) + RotatedAxis()

FeaturePlot(object = livMets, features = c("percent.mt"))
FeaturePlot(object = livMets, features = c("percent.ribo"))
FeaturePlot(object = livMets, features = c("nCount_RNA"))
FeaturePlot(object = livMets, features = c("nFeature_RNA"))


# selecting cell groups to assign clusters manually
DimPlot(livMets, reduction = "umap",label = T,pt.size = 2)
# https://satijalab.org/seurat/articles/visualization_vignette.html
plot <- DimPlot(livMets, reduction = "umap")
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(livMets, cells = select.cells) <- 24
# Idents(livMets)
DimPlot(livMets,label = T, pt.size = 2)
table((livMets@active.ident))
levels(livMets@active.ident)
livMets@active.ident <- factor(livMets@active.ident, levels = c("0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19","20","21","22","23","24"))
DimPlot(livMets,label = T, pt.size = 1)
DotPlot(livMets, features = markerGenes) + RotatedAxis()

# get cluster markers
markers <- FindAllMarkers(livMets,only.pos = T, logfc.threshold = 0.25)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
head(markers,20)
# write.csv(markers, file = "code/Collagen_DDR1_PDACmets/results/markers_pdac_LiM_clusters.csv")

# giving cell types
DimPlot(livMets,label = T)
# livMets[["Cell.clusters"]] <- Idents(object = livMets)
DimPlot(livMets,label = T,group.by = "Cell.clusters")
livMets <- SetIdent(livMets, value = "Cell.clusters")
livMets <- RenameIdents(object = livMets, 
                           '2' = "Epithelial-Met",
                           '6' = "Epithelial-Met",
                           '10' = "Epithelial-Met",
                           '8' = "Epithelial-Met-cycling",
                           
                           '17' = "Cholangiocytes",
                           '16' = "Hepatocytes1",
                           '20' = "Hepatocytes2",
                           
                           '13' = "VSMC",
                           '21' = "Fibroblasts",
                           '7' = "Endothelial",
                           
                           '5' = "Myeloid",
                           '9' = "Myeloid",
                           '11' = "Myeloid",
                           '23' = "Myeloid-cycling",
                        
                           '19' = "DC",
                           
                           '3' = "B",
                           '14' = "Plasma",
                           '22' = "Plasma",
                           
                           '0' = "T_NK",
                           '1' = "T_NK",
                           '4' = "T_NK",
                           '12' = "T_NK",
                           '15' = "T_NK",
                           '24' = "T_NK",
                           
                           '18' = "Cycling")


livMets[["Cell.types"]] <- Idents(object = livMets)
DimPlot(livMets,label = T, group.by = "Cell.types")
table(livMets@meta.data$Cell.types)
table(livMets@meta.data$Cell.clusters)
table(livMets@meta.data$cell.type)

DimPlot(livMets, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(livMets, reduction = "umap",group.by = "Cell.clusters", label = F)

DotPlot(livMets, features = markerGenes) + RotatedAxis()


# giving cell types- simpler
DimPlot(livMets,label = T)
# livMets[["Cell.clusters"]] <- Idents(object = livMets)
DimPlot(livMets,label = T,group.by = "Cell.clusters")
livMets <- SetIdent(livMets, value = "Cell.clusters")
livMets <- RenameIdents(object = livMets, 
                        '2' = "Epithelial-Met",
                        '6' = "Epithelial-Met",
                        '10' = "Epithelial-Met",
                        '8' = "Epithelial-Met",
                        
                        '17' = "Cholangiocytes",
                        '16' = "Hepatocytes",
                        '20' = "Hepatocytes",
                        
                        '13' = "VSMC",
                        '21' = "Fibroblasts",
                        '7' = "Endothelial",
                        
                        '5' = "Myeloid",
                        '9' = "Myeloid",
                        '11' = "Myeloid",
                        '23' = "Myeloid",
                        # '22' = "Myeloid",
                        '19' = "DC",
                        
                        '3' = "B",
                        '14' = "Plasma",
                        '22' = "Plasma",
                        
                        '0' = "T_NK",
                        '1' = "T_NK",
                        '4' = "T_NK",
                        '12' = "T_NK",
                        '15' = "T_NK",
                        '24' = "T_NK",
                        
                        # '5' = "NK",
                        '18' = "Cycling")


livMets[["Cell.types.simple"]] <- Idents(object = livMets)
DimPlot(livMets,label = T, group.by = "Cell.types.simple")
table(livMets@meta.data$Cell.types)
table(livMets@meta.data$Cell.clusters)
table(livMets@meta.data$cell.type)

DimPlot(livMets, reduction = "umap",group.by = "Cell.types.simple", label = F)
DimPlot(livMets, reduction = "umap",group.by = "Cell.clusters", label = F)
# saveRDS(livMets, file = "results/seuratObj_pdac_LiM.rds")


# **************************************** M1/M2 analysis
# M1/M2 signatures were obtained from ref: PMIDs: 35385733, 30721157
m1Sig <- toupper(unlist(strsplit("Azin1, Cd38, Cxcl10, Cxcl9, Fpr2, Il18, Il1b, Irf5, Nifkbiz, Tlr4, Tnf, Cd80",split=", ")))
m2Sig <- toupper(unlist(strsplit("Alox5, Arg1, Chil3, Cd163, Il10, Il10ra, Il10rb, Irf4, Kif4, Mrc1, Myc, Socs2, Tgm2",split=", ")))

cells = WhichCells(object = livMets, idents = "Myeloid")

genes <- m1Sig
genes <- intersect(genes,rownames(livMets@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = livMets, features = c(gene), max.cutoff = "q99",pt.size = 1, order = T)


DimPlot(livMets, reduction = "umap",group.by = "Cell.clusters", label = T,cells = cells)

sig <- m1Sig
sigName <- "m1Sig"
callName <- paste0(sigName,"1")
livMets <- AddModuleScore(livMets,features = list(sig), name = sigName)
FeaturePlot(livMets, callName, max.cutoff = "q95",pt.size = 1,order = T,cells = cells) + scale_color_viridis_c(option = "viridis") +labs(title = "Human PDAC LiM M1 signature")
VlnPlot(livMets, group.by = "Cell.clusters", features = callName, idents = "Myeloid")

genes <- m2Sig
genes <- intersect(genes,rownames(livMets@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = livMets, features = c(gene), max.cutoff = "q99",pt.size = 1, order = T)

sig <- m2Sig
sigName <- "m2Sig"
callName <- paste0(sigName,"1")
livMets <- AddModuleScore(livMets,features = list(sig), name = sigName)
FeaturePlot(livMets, callName, max.cutoff = "q99",pt.size = 1,order = T,cells = cells) + scale_color_viridis_c(option = "viridis") +labs(title = "Human PDAC LiM M2 signature")
VlnPlot(livMets, group.by = "Cell.clusters", features = callName, idents = "Myeloid")

DimPlot(livMets,label = T,group.by = "Cell.clusters")
livMets <- SetIdent(livMets, value = "Cell.clusters")
livMets <- RenameIdents(object = livMets, 
                        '2' = "Epithelial-Met",
                        '6' = "Epithelial-Met",
                        '10' = "Epithelial-Met",
                        '8' = "Epithelial-Met",
                        
                        '17' = "Cholangiocytes",
                        '16' = "Hepatocytes",
                        '20' = "Hepatocytes",
                        
                        '13' = "VSMC",
                        '21' = "Fibroblasts",
                        '7' = "Endothelial",
                        
                        '5' = "M1-like",
                        '9' = "M1-like",
                        '11' = "M2-like",
                        '23' = "M1-like",
                        '19' = "DC",
                        
                        '3' = "B",
                        '14' = "Plasma",
                        '22' = "Plasma",
                        
                        '0' = "T_NK",
                        '1' = "T_NK",
                        '4' = "T_NK",
                        '12' = "T_NK",
                        '15' = "T_NK",
                        '24' = "T_NK",
                        
                        '18' = "Cycling")



livMets[["Cell.types.simple.M1M2"]] <- Idents(object = livMets)
DimPlot(livMets,label = T, group.by = "Cell.types.simple.M1M2")
table(livMets@meta.data$Cell.types)
table(livMets@meta.data$Cell.clusters)
table(livMets@meta.data$cell.type)

DimPlot(livMets, reduction = "umap",group.by = "Cell.types.simple", label = F)
DimPlot(livMets, reduction = "umap",group.by = "Cell.types.simple.M1M2", label = F)
DimPlot(livMets, reduction = "umap",group.by = "Cell.clusters", label = F)
# saveRDS(livMets, file = "results/seuratObj_pdac_LiM.rds")


# getting the cell type annotation
colnames(livMets@meta.data)
head(livMets@meta.data)
celltypes <- livMets@meta.data[,c(2,3,4,7,8,9,10,13)]
head(celltypes)
# write.csv(celltypes, file = "results/cellAnnotation_pdac_LiM_wliver.csv")

# **************************************** working with saved  data
livMets <- readRDS(file = "results/seuratObj_pdac_LiM.rds")
DimPlot(livMets, label = T)

DimPlot(livMets, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(livMets, reduction = "umap",group.by = "Cell.types.simple", label = T)
DimPlot(livMets, reduction = "umap",group.by = "Cell.types.simple.M1M2", label = F)
DimPlot(livMets, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(livMets, reduction = "umap",group.by = "orig.ident", label = T)


#******************************************************#******************************************************
#******************************************************PLOTS 
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")


show_col(pal_ucscgb("default")(20))
pal_ucscgb("default")(20)

livMets <- SetIdent(livMets, value = "Cell.types")
req_subset <- livMets
levels(req_subset@active.ident)


req_subset <- SetIdent(req_subset, value = "Cell.types.simple")
# req_subset <- livMets
levels(req_subset@active.ident)

cols <- c("#FF0000FF","#9900CCFF","#00468BFF","#FFCC00FF","#1B1919FF","#CCFF00FF","#FDBF6F","#58593FFF","#42B540FF","#B2df8a","#AD002AFF","#CC99FFFF")
DimPlot(object = req_subset, group.by = "Cell.types.simple", label = T, cols = cols)

# DimPlot(object = req_subset, group.by = "Cell.types.simple", label = T) +   scale_color_ucscgb()  + labs(title = "Human PDAC LiM")
DimPlot(object = req_subset, group.by = "Cell.types.simple", label = T, cols = cols) + labs(title = "Human PDAC LiM")
ggsave(filename = paste0("results/figs/hum_PDACLiM_celltypes_label.pdf"))

# DimPlot(object = req_subset, group.by = "Cell.types.simple", label = F) +   scale_color_ucscgb()  + labs(title = "Human PDAC LiM")
DimPlot(object = req_subset, group.by = "Cell.types.simple", label = F, cols = cols) + NoLegend()  + labs(title = "Human PDAC LiM")
ggsave(filename = paste0("results/figs/hum_PDACLiM_celltypes_nolegend.pdf"))


# gene umaps
DefaultAssay(req_subset) <- "RNA"
# dOTPLOTS
markers <- c("EPCAM","KRT8",# "
             "ANXA4",'SOX9',
             "ALB","APOE",#"APOB","CPS1",
             "MYH11","ACTA2",#,"ACTG2"
             "COLEC11","DCN","COL1A1",#"COL3A1","LOX","TIMP1",
             "KDR","PECAM1",
             "CD68","FCGR3A", "LYZ",
             "FCER1A","CPA3",# 
             "MS4A1","CD79A",#"IGKC",
             "JCHAIN",
             "IL7R","CD3D",#,"CD8A""TRAC","FOXP3",
             "NKG7","GNLY",
             "TOP2A","MKI67")#,#"CENPF",)
             
DotPlot(req_subset, features = markers) + RotatedAxis() & scale_color_distiller(palette = "RdBu")
ggsave(filename = "results/figs/hum_PDACLiM_dotplot_markers.pdf", width = 11, height = 6)

genes <- unlist(strsplit("EPCAM, KRT19, DDR1, DDR2, LAIR1, MRC2, ITGB1, NFE2L2, SDC1, TFAM, CDC42", split = ", "))
genes <- intersect(genes,rownames(req_subset@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q99",pt.size = 1, order = T)

for (gene in genes) {
  FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q95",pt.size = 1, order = F) #+ scale_color_viridis_c()
  ggsave(filename = paste0("results/figs/hum_PDACLiM_gene_",gene,".pdf"))
}


genes <- unlist(strsplit("MMP1,MMP2,MMP8,MMP9,MMP13,MMP14", split = ","))
genes <- intersect(genes,rownames(req_subset@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q99",pt.size = 1, order = F)

for (gene in genes) {
  FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q95",pt.size = 1, order = F) #+ scale_color_viridis_c()
  ggsave(filename = paste0("results/figs/hum_PDACLiM_gene_",gene,".pdf"))
}




# **************************************** M1/M2
req_subset <- SetIdent(req_subset, value = "Cell.types.simple") 
cells = WhichCells(object = req_subset, idents = "Myeloid")


DimPlot(object = req_subset, group.by = "Cell.types.simple", label = F, cols = cols[7], cells = cells, pt.size = 3)+ xlim(-11,-4) +ylim(-12,-5)
ggsave(filename = paste0("results/figs/hum_PDACLiM_M1M2_celltypes.pdf"))

req_subset <- SetIdent(req_subset, value = "Cell.types.simple.M1M2")
col1 <- c("#FF9900FF", "#6699FFFF")
DimPlot(object = req_subset, group.by = "Cell.types.simple.M1M2", label = F, cols = col1, cells = cells, pt.size = 3)+ xlim(-11,-4) +ylim(-12,-5)
ggsave(filename = paste0("results/figs/hum_PDACLiM_M1M2_celltypes_M1M2.pdf"))


req_subset <- SetIdent(req_subset, value = "Cell.types.simple") 
sigName <- "m1Sig"
callName <- paste0(sigName,"1")
VlnPlot(req_subset, group.by = "Cell.types.simple.M1M2", features = callName, idents = "Myeloid",cols = col1)
ggsave(filename = paste0("results/figs/hum_PDACLiM_M1M2_M1sig.pdf"), width = 5, height = 8)

sigName <- "m2Sig"
callName <- paste0(sigName,"1")
VlnPlot(req_subset, group.by = "Cell.types.simple.M1M2", features = callName, idents = "Myeloid",cols = col1)
ggsave(filename = paste0("results/figs/hum_PDACLiM_M1M2_M2sig.pdf"), width = 5, height = 8)


genes <- unlist(strsplit("MMP1,MMP2,MMP8,MMP9,MMP13,MMP14", split = ","))
genes <- intersect(genes,rownames(req_subset@assays$RNA@data))
gene <- genes[1]
FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q99", order = T, cells = cells, pt.size = 3)+ xlim(-11,-4) +ylim(-12,-5)

for (gene in genes) {
  FeaturePlot(object = req_subset, features = c(gene), max.cutoff = "q95", order = T, cells = cells, pt.size = 3)+ xlim(-11,-4) +ylim(-12,-5)
  ggsave(filename = paste0("results/figs/hum_PDACLiM_M1M2_gene_",gene,".pdf"))
}
