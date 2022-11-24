library(Seurat);library("RColorBrewer")

pbmc.data=read.csv("matrix/E07.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE07.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E07.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E07.5")
pbmc.data=read.csv("matrix/E08.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE08.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E08.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E08.5")
pbmc.data=read.csv("matrix/E09.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE09.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E09.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E09.5")
pbmc.data=read.csv("matrix/E10.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE10.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E10.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E10.5")
pbmc.data=read.csv("matrix/E10.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE10.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
E10.5h <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E10.5h")
pbmc.data=read.csv("matrix/E11.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE11.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
E11.5h <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E11.5h")
pbmc.data=read.csv("matrix/E12.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE12.5h.csv",check.names = FALSE, header=TRUE,row.names=1)
E12.5h <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E12.5h")
pbmc.data=read.csv("matrix/E12.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE12.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E12.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E12.5")
pbmc.data=read.csv("matrix/E13.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE13.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E13.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E13.5")
pbmc.data=read.csv("matrix/E14.5.csv",check.names = FALSE, header=TRUE,row.names=1)
meta.data=read.csv("matrix/dE14.5.csv",check.names = FALSE, header=TRUE,row.names=1)
E14.5 <- CreateSeuratObject(counts = pbmc.data,meta.data =meta.data, project = "E14.5")

pbmc <- merge(x = E07.5, y = c(E08.5,E09.5,E10.5, E10.5h, E11.5h, E12.5h,  E12.5, E13.5,E14.5 ),  
              add.cell.ids = c("E07.5","E08.5","E09.5","E10.5", "E10.5h", "E11.5h", "E12.5h","E12.5",  "E13.5","E14.5"), project = "mpl")
table(Idents(pbmc))

pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data =pbmc@meta.data, min.cells = 3, project = "mtr")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(pbmc))
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,verbose = F)
DimPlot(pbmc, reduction = "pca",label=T)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:10,umap.method ='umap-learn',metric= 'correlation')
DimPlot(pbmc,reduction = "umap",  label=T) 

pbmc <- FindNeighbors(object = pbmc, dims = 1:13)
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.05)
table(Idents(pbmc))
DimPlot(pbmc, label=T) 



#plot <- DimPlot(pbmc, reduction = "umap")
#s1 <- CellSelector(plot = plot)
#s2 <- CellSelector(plot = plot)
#write.table(s1, file="cluster_E_cell.txt",row.names = F, col.names = F,quote = F)  
#write.table(s2, file="cluster_F_cell.txt",row.names = F, col.names = F,quote = F)  
s1=readLines("cluster_E_cell.txt")
s2=readLines("cluster_F_cell.txt")


Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "N"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "M"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "L"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "K"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "J"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "I"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "I"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "H"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "G"
Idents(pbmc, cells = s2) <- "F"
Idents(pbmc, cells = s1) <- "E"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "D"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "C"
Idents(pbmc, cells = WhichCells(pbmc,idents = "0")) <- "B"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "A"
pbmc$cluster_14=Idents(pbmc)
DimPlot(object = pbmc,label = T,cols=colorRampPalette(brewer.pal(8, "Set1"))(14))

new.cluster.ids <- c("SpA-TGC ","Trophoblast cell","P-TGC","Primitive endoderm cell","Yolk sac epithelial cell","Embryo stem cell","Decidual stromal cell",
                     "Decidual pericyte","Endothelial cell","Embryo stromal cell","Immune cell","Megakaryocyte","Erythrocyte","Hematopoietic cell")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$cell_type_14=Idents(pbmc)
DimPlot(object = pbmc,label = T,cols=colorRampPalette(brewer.pal(8, "Set1"))(14))

pbmc$orig=pbmc$stage
pbmc$sample=pbmc$orig.ident

Idents(pbmc)='orig.ident'
table(Idents(pbmc))
new.cluster.ids <- c("E07.5","E08.5","E09.5","E10.5", "E10.5", "E11.5", "E12.5","E12.5",  "E13.5","E14.5")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stage=Idents(pbmc)
DimPlot(object = pbmc,label = T)


Idents(pbmc)='orig.ident'
table(Idents(pbmc))
new.cluster.ids <- c("E07.5-E08.5","E07.5-E08.5","E09.5-E10.5","E09.5-E10.5", "E09.5-E10.5", "E11.5-E14.5", "E11.5-E14.5","E11.5-E14.5",  "E11.5-E14.5","E11.5-E14.5")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stages=Idents(pbmc)
DimPlot(object = pbmc,label = T)

Idents(pbmc)="cluster_14"  
DimPlot(object = pbmc,label = T)

saveRDS(pbmc, file = "MPL39603.rds")
write.csv(pbmc@meta.data, file="MPL39603-meta.csv")  




pbmc=readRDS("MPL39603.rds");pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
Idents(pbmc)="cluster_14"
pbmcx=subset(pbmc,idents='B');pbmcx
pbmcx <- subset(x = pbmcx, subset =nFeature_RNA <= 2000 | nCount_RNA >= 60000  | percent.mt<=1);pbmcx
low=colnames(pbmcx)

#plot <- DimPlot(pbmc, reduction = "umap")
#extend <- CellSelector(plot = plot)
extend=readLines("Multiplet-extend.txt")

pbmcx=subset(pbmc,idents=c('A','B','C'));pbmcx
Idents(pbmcx)='doublet'
table(Idents(pbmcx))
scru=WhichCells(pbmcx,idents = "1")

Idents(pbmcx)="cluster_14"
Idents(pbmcx, cells = c(low,extend,scru)) <- "Multiple and Low quality cells"
table(Idents(pbmcx))
DimPlot(object = pbmcx,label = T)


pbmc <- SubsetData(pbmcx, ident.use = c('A','B','C'));pbmc
pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data =pbmc@meta.data, min.cells = 3, project = "mtr");pbmc
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(pbmc))
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,verbose = F)
DimPlot(pbmc, reduction = "pca", label=T)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:11,umap.method ='umap-learn',metric= 'correlation')
DimPlot(pbmc, label=T) 



pbmc <- FindNeighbors(object = pbmc, dims = 1:11)
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.6)
pbmc <- ReorderIdent(pbmc, c("PC_5"), reorder.numeric = T,reverse=T)
new.cluster.ids <- c(1,2,17,3,6,4,7,9,8,13,11,10,14,5, 12,15,16)
names(x = new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
DimPlot(pbmc, label=T) 
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14"
Idents(pbmc, cells = WhichCells(pbmc,idents = "13")) <- "13"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "7"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1"
pbmc$cluster_17=Idents(pbmc)
DimPlot(pbmc, label=T) 


C18=readLines("cluster_18_cell.txt")
C19=readLines("cluster_19_cell.txt")
Idents(pbmc, cells = C19) <- "19"
Idents(pbmc, cells = C18) <- "18"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14"
Idents(pbmc, cells = WhichCells(pbmc,idents = "13")) <- "13"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "7"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1"
pbmc$cluster_19=Idents(pbmc)
C19=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1:9,22:24,18,11:14,27,30)]
DimPlot(pbmc, label=T,cols=C19,group.by='cluster_19') 

Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18 EPC Migratory Cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "15 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "13")) <- "13 EPC cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12 Secondary P-TGC \n\     Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11 Secondary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10 Primary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9   S-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "8   S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "7   S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6   S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "5   SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4   SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3   LaTP 2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2   LaTP"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1   TSC and ExE cell"
pbmc$type_19=Idents(pbmc)
DimPlot(pbmc, label=T,cols=C19) 

Idents(pbmc)="cluster_19"  
DimPlot(pbmc,label = T,cols=C19)


pbmc2=pbmc

Idents(pbmc2) <- "cluster_19"
pbmc <- subset(pbmc2, idents = c(6,7,14,15));pbmc
DimPlot(pbmc, label=T)
pbmc <- subset(pbmc2, idents = c(6,7));pbmc
DimPlot(pbmc, label=T)

pbmc <- FindNeighbors(pbmc,  dims = 1:9)
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.4) 
table(Idents(pbmc));DimPlot(object = pbmc,label = T)
Idents(pbmc2, cells = WhichCells(pbmc,idents = "2")) <- "P3"
Idents(pbmc2, cells = WhichCells(pbmc,idents = "4")) <- "P2"

pbmc <- subset(pbmc2, idents = c(13));pbmc
DimPlot(pbmc, label=T)
pbmc <- FindNeighbors(pbmc,  dims = 1:11)
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.2) 
table(Idents(pbmc));DimPlot(object = pbmc,label = T )

Idents(pbmc2, cells = WhichCells(pbmc,idents = "1")) <- "P1"
Idents(pbmc2, cells = WhichCells(pbmc,idents = "2")) <- "E2"
Idents(pbmc2, cells = WhichCells(pbmc,idents = "0")) <- "E1"


DimPlot(object = pbmc2,label = T )
pbmc=pbmc2

Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "7"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1"
pbmc$cluster_17E=Idents(pbmc)
DimPlot(object = pbmc,label = T )

table(Idents(pbmc))


Idents(pbmc, cells = WhichCells(pbmc,idents = "P3")) <- "P3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P2")) <- "P2 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1")) <- "P1 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2 The progenitor of S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1 The progenitor of SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18 EPC Migratory Cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "15 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12 Secondary P-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11 Secondary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10 Primary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9   S-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "8   S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "7   S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6   S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "5   SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4   SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3   LaTP 2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2   LaTP"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1   TSC and ExE cell"
pbmc$type_17E=Idents(pbmc)
DimPlot(object = pbmc,label = T )

Idents(pbmc)='cluster_17E'
Idents(pbmc, cells = WhichCells(pbmc,idents = "P3")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P2")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "6-7-8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "6-7-8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6-7-8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "4-5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4-5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1"
pbmc$cluster_21=Idents(pbmc)
DimPlot(object = pbmc,label = T )

Idents(pbmc)='cluster_17E'
Idents(pbmc, cells = WhichCells(pbmc,idents = "P3")) <- "P1-2-3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P2")) <- "P1-2-3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1")) <- "P1-2-3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2 The progenitor of S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1 The progenitor of SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18 EPC Migratory Cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "14-15 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14-15 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12 Secondary P-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11 Secondary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10 Primary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9   S-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "8")) <- "6-7-8 S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "7")) <- "6-7-8 S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6")) <- "6-7-8 S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "4-5 SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4")) <- "4-5 SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3   LaTP 2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2   LaTP"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1   TSC and ExE cell"

pbmc$type_21=Idents(pbmc)
DimPlot(object = pbmc,label = T )

Idents(pbmc)='cluster_17E'
new.cluster.ids <- c("TSC and ExE cell","Chorion branch","Chorion branch","Chorion branch","Chorion branch",
                     "Sinusoid branch","Sinusoid branch","Sinusoid branch","Sinusoid branch",
                     "P-TGC branch","P-TGC branch","P-TGC branch", 
                     "Spongio-branch", "Spongio-branch","Spongio-branch","Spongio-branch",
                     "P-TGC branch","Chorion branch","Spongio-branch","Sinusoid branch",
                     "Biopotential progenitor","Biopotential progenitor","Biopotential progenitor")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
Idents(pbmc, cells = WhichCells(pbmc,idents = "P-TGC branch")) <- "P-TGC branch"
Idents(pbmc, cells = WhichCells(pbmc,idents = "Biopotential progenitor")) <- "Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "TSC and ExE cell")) <- "TSC and ExE cell"
table(Idents(pbmc))
pbmc$branch=Idents(pbmc)
DimPlot(object = pbmc,label = T,cols=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1,18,23,4,8,13)] )

cc.genes <- readLines(con = "regev_lab_cell_cycle_gene.txt");   s.genes <- cc.genes[1:43];    g2m.genes <- cc.genes[44:97]
pbmc <- CellCycleScoring(object = pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(pbmc)

Idents(pbmc)='cluster_17E'
DimPlot(object = pbmc,label = T )
saveRDS(pbmc,"MTR15682.rds");pbmc
write.csv(pbmc@meta.data, file="MTR15682-meta.csv")  








sessionInfo()

R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936    LC_MONETARY=Chinese (Simplified)_China.936
[4] LC_NUMERIC=C                               LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] RColorBrewer_1.1-2 Seurat_3.1.5      

loaded via a namespace (and not attached):
[1] tsne_0.1-3          nlme_3.1-139        bitops_1.0-6        fs_1.3.1            usethis_1.5.0       devtools_2.0.2     
[7] RcppAnnoy_0.0.12    httr_1.4.0          rprojroot_2.0.2     sctransform_0.2.0   tools_3.6.0         R6_2.4.0           
[13] irlba_2.3.3         KernSmooth_2.23-15  uwot_0.1.8          lazyeval_0.2.2      colorspace_1.4-1    npsurv_0.4-0       
[19] withr_2.1.2         gridExtra_2.3       tidyselect_1.1.0    prettyunits_1.0.2   processx_3.3.1      compiler_3.6.0     
[25] cli_1.1.0           desc_1.2.0          plotly_4.9.0        labeling_0.3        caTools_1.17.1.2    scales_1.0.0       
[31] lmtest_0.9-37       ggridges_0.5.1      callr_3.2.0         pbapply_1.4-0       stringr_1.4.0       digest_0.6.18      
[37] pkgconfig_2.0.2     htmltools_0.3.6     sessioninfo_1.1.1   htmlwidgets_1.3     rlang_0.4.7         rstudioapi_0.10    
[43] generics_0.0.2      zoo_1.8-7           jsonlite_1.6        ica_1.0-2           gtools_3.8.1        dplyr_1.0.2        
[49] magrittr_1.5        patchwork_1.0.0     Matrix_1.2-18       Rcpp_1.0.5          munsell_0.5.0       ape_5.3            
[55] reticulate_1.12     lifecycle_0.2.0     stringi_1.4.3       MASS_7.3-51.4       pkgbuild_1.0.3      gplots_3.0.1.1     
[61] Rtsne_0.15          plyr_1.8.4          grid_3.6.0          parallel_3.6.0      gdata_2.18.0        listenv_0.7.0      
[67] ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-38     cowplot_0.9.4       splines_3.6.0       ps_1.3.0           
[73] pillar_1.4.6        igraph_1.2.6        reshape2_1.4.3      future.apply_1.2.0  codetools_0.2-16    pkgload_1.0.2      
[79] leiden_0.3.5        glue_1.4.2          lsei_1.2-0          data.table_1.12.2   remotes_2.0.4       BiocManager_1.30.4 
[85] png_0.1-7           vctrs_0.3.4         testthat_2.3.2      gtable_0.3.0        RANN_2.6.1          purrr_0.3.2        
[91] tidyr_1.1.2         future_1.13.0       assertthat_0.2.1    ggplot2_3.1.1       rsvd_1.0.0          viridisLite_0.3.0  
[97] survival_2.44-1.1   tibble_3.0.3        memoise_1.1.0       cluster_2.0.8       globals_0.12.4      fitdistrplus_1.0-14
[103] ellipsis_0.3.0      ROCR_1.0-7



