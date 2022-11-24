library(Seurat)
PBMC <- Read10X(data.dir = "filtered_feature_bc_matrix");dim(PBMC)
pbmc <- CreateSeuratObject(counts = PBMC,  project = "MPL",min.cells = 3,names.field = 2, names.delim = "-" )

Idents(pbmc)="orig.ident"
new.cluster.ids <- c('E7.5','E8.5','E9.5','E11.5','E13.5')
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stage=Idents(pbmc)
table(Idents(pbmc))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(object = pbmc,pt.size=0.1, group.by = "stage",features = c("nFeature_RNA","percent.mt")) 

pbmc <- subset(pbmc, subset = nFeature_RNA >1000 & percent.mt < 10);pbmc
table(Idents(pbmc))
pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data, min.cells = 3, project = "mtr",names.field = 2, names.delim = "-");pbmc

write.csv(as.matrix(subset(pbmc,idents =1)$RNA@counts), "V3E7.csv"); write.csv(subset(pbmc,idents =1)@meta.data, "V3E7D.csv")
write.csv(as.matrix(subset(pbmc,idents =2)$RNA@counts), "V3E8.csv"); write.csv(subset(pbmc,idents =2)@meta.data, "V3E8D.csv")
write.csv(as.matrix(subset(pbmc,idents =3)$RNA@counts), "V3E9.csv"); write.csv(subset(pbmc,idents =3)@meta.data, "V3E9D.csv")
write.csv(as.matrix(subset(pbmc,idents =4)$RNA@counts), "V3E11.csv"); write.csv(subset(pbmc,idents =4)@meta.data, "V3E11D.csv")
write.csv(as.matrix(subset(pbmc,idents =5)$RNA@counts), "V3E13.csv"); write.csv(subset(pbmc,idents =5)@meta.data, "V3E13D.csv")


library(Seurat)
library("RColorBrewer")
M1=read.csv("V3E7.csv",check.names = FALSE, header=TRUE,row.names=1)
M2=read.csv("V3E8.csv",check.names = FALSE, header=TRUE,row.names=1)
M3=read.csv("V3E9.csv",check.names = FALSE, header=TRUE,row.names=1)
M4=read.csv("V3E11.csv",check.names = FALSE, header=TRUE,row.names=1)
M5=read.csv("V3E13.csv",check.names = FALSE, header=TRUE,row.names=1)
meta=read.csv("V3E7D.csv",check.names = FALSE, header=TRUE,row.names=1);  M1<- CreateSeuratObject(M1, meta.data =meta, project = "V3E7")
meta=read.csv("V3E8D.csv",check.names = FALSE, header=TRUE,row.names=1);  M2<- CreateSeuratObject(M2, meta.data =meta, project = "V3E8")
meta=read.csv("V3E9D.csv",check.names = FALSE, header=TRUE,row.names=1);  M3<- CreateSeuratObject(M3, meta.data =meta, project = "V3E9")
meta=read.csv("V3E11D.csv",check.names = FALSE, header=TRUE,row.names=1); M4<- CreateSeuratObject(M4, meta.data =meta, project = "V3E11")
meta=read.csv("V3E13D.csv",check.names = FALSE, header=TRUE,row.names=1); M5<- CreateSeuratObject(M5, meta.data =meta, project = "V3E13")

pbmc <- merge(x =M1, y = c(M2,M3,M4,M5), 
              add.cell.ids = c("E7.5-V3","E8.5-V3","E9.5-V3","E11.5-V3","E13.5-V3"), project = "mtr");      pbmc

pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data, min.cells = 3, project = "MF");pbmc
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(pbmc))
pbmc <- ScaleData(pbmc); pbmc <- RunPCA(pbmc,verbose = F)

Idents(pbmc)='doublet'
pbmc=subset(pbmc,idents='0')
MF=pbmc;pbmc

pbmcx <- readRDS("MPL39603.rds")
Idents(pbmcx)='doublet'
pbmcx=subset(pbmcx,idents='0');pbmcx
MM <- RunUMAP(pbmcx, reduction = "pca", dims = 1:10,verbose=F, return.model = TRUE)
Idents(MM)="cluster_14"
DimPlot(MM, label=T) 
MM$umap@cell.embeddings=pbmcx$umap@cell.embeddings
MM$umap@misc$model$embedding=pbmcx$umap@cell.embeddings
DimPlot(MM, label=T) 

anchors <- FindTransferAnchors(
  reference= MM,
  query = MF,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
)
MF <- MapQuery(
  anchorset = anchors,
  query = MF,
  reference = MM,
  refdata = list(
    celltype1 = "cluster_14"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)
DimPlot(MF, label = T,reduction = "ref.umap", group.by = "predicted.celltype1")

Idents(MF)="predicted.celltype1"
Idents(MF, cells = WhichCells(MF,idents = "N")) <- "N"
Idents(MF, cells = WhichCells(MF,idents = "M")) <- "M"
Idents(MF, cells = WhichCells(MF,idents = "L")) <- "L"
Idents(MF, cells = WhichCells(MF,idents = "K")) <- "K"
Idents(MF, cells = WhichCells(MF,idents = "J")) <- "J"
Idents(MF, cells = WhichCells(MF,idents = "I")) <- "I"
Idents(MF, cells = WhichCells(MF,idents = "H")) <- "H"
Idents(MF, cells = WhichCells(MF,idents = "G")) <- "G"
Idents(MF, cells = WhichCells(MF,idents = "F")) <- "F"
Idents(MF, cells = WhichCells(MF,idents = "E")) <- "E"
Idents(MF, cells = WhichCells(MF,idents = "D")) <- "D"
Idents(MF, cells = WhichCells(MF,idents = "C")) <- "C"
Idents(MF, cells = WhichCells(MF,idents = "B")) <- "B"
Idents(MF, cells = WhichCells(MF,idents = "A")) <- "A"
MF$cluster_14=Idents(MF)
pbmc=MF
DimPlot(pbmc, label = T,reduction = "ref.umap",cols=colorRampPalette(brewer.pal(8, "Set1"))(14))
new.cluster.ids <- c("SpA-TGC ","Trophoblast cell","P-TGC","Primitive endoderm cell","Yolk sac epithelial cell","Embryo stem cell","Decidual stromal cell",
                     "Decidual pericyte","Endothelial cell","Embryo stromal cell","Immune cell","Megakaryocyte","Erythrocyte","Hematopoietic cell")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$type_14=Idents(pbmc)
DimPlot(pbmc, label = T,reduction = "ref.umap",cols=colorRampPalette(brewer.pal(8, "Set1"))(14))




Idents(pbmc)="cluster_14"
pbmcy=subset(pbmc,idents="B");pbmcy
pbmcy <- subset(pbmcy, subset = nFeature_RNA <=2000 | percent.mt <= 1);pbmcy
Idents(pbmc, cells = WhichCells(pbmcy)) <- "Low"
MF=subset(pbmc,idents=c("A","B","C"));MF

MF <- CreateSeuratObject(MF$RNA@counts, meta.data = MF@meta.data, min.cells = 3, project = "MF");MF
MF <- NormalizeData(MF)
MF <- FindVariableFeatures(MF, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(MF))
MF <- ScaleData(MF); MF <- RunPCA(MF,verbose = F)

pbmcx <- readRDS("MTR15682.rds")
MM <- RunUMAP(pbmcx, reduction = "pca", dims = 1:11,verbose=F, return.model = TRUE)
Idents(MM)="cluster_19"
MM$umap@cell.embeddings=pbmcx$umap@cell.embeddings
MM$umap@misc$model$embedding=pbmcx$umap@cell.embeddings
DimPlot(MM, label=T)

anchors <- FindTransferAnchors(
  reference= MM,
  query = MF,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
)
MF <- MapQuery(
  anchorset = anchors,
  query = MF,
  reference = MM,
  refdata = list(
    C19 = "cluster_19",
    C17E = "cluster_17E",
    C21 = "cluster_21",
    type = "type_21"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)
DimPlot(MF, label = T,reduction = "ref.umap", group.by = "predicted.celltype1")

pbmc=MF
Idents(pbmc)="predicted.C19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18"
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
Idents(pbmc)="predicted.C17E"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P3")) <- "P3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P2")) <- "P2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1")) <- "P1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1"
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

Idents(pbmc)="predicted.C21"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1-2-3")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14-15")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6-7-8")) <- "6-7-8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4-5")) <- "4-5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1"
pbmc$cluster_21=Idents(pbmc)
C21=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1:4,7,9,22:24,12:14,27,30,15,10,18)]
DimPlot(pbmc,label=T, cols=C21,reduction='ref.umap')
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1-2-3")) <- "P1-2-3 EPC progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2 The progenitor of S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1 The progenitor of SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18 EPC Migratory Cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14-15")) <- "14-15 SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12 Secondary P-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11 Secondary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10 Primary P-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "9")) <- "9   S-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "6-7-8")) <- "6-7-8 S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "4-5")) <- "4-5 SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "3")) <- "3   LaTP 2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "2")) <- "2   LaTP"
Idents(pbmc, cells = WhichCells(pbmc,idents = "1")) <- "1   TSC and ExE cell"
pbmc$type_21=Idents(pbmc)
DimPlot(pbmc,label=T, cols=C21,reduction='ref.umap')







sessionInfo()

R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /home/jxx/anaconda3/envs/r4/lib/libmkl_rt.so.1

locale:
[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
[9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] RColorBrewer_1.1-2 SeuratObject_4.0.0 Seurat_4.0.0      

loaded via a namespace (and not attached):
[1] nlme_3.1-151          matrixStats_0.57.0    RcppAnnoy_0.0.18     
[4] httr_1.4.2            repr_1.1.0            sctransform_0.3.2    
[7] tools_4.0.3           utf8_1.1.4            R6_2.5.0             
[10] irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-18   
[13] uwot_0.1.10           mgcv_1.8-33           DBI_1.1.0            
[16] lazyeval_0.2.2        colorspace_2.0-0      tidyselect_1.1.2     
[19] gridExtra_2.3         compiler_4.0.3        cli_3.2.0            
[22] Cairo_1.5-12.2        plotly_4.9.2.2        labeling_0.4.2       
[25] scales_1.1.1          spatstat.data_1.7-0   lmtest_0.9-38        
[28] ggridges_0.5.2        pbapply_1.4-3         goftest_1.2-2        
[31] spatstat_1.64-1       pbdZMQ_0.3-4          stringr_1.4.0        
[34] digest_0.6.27         spatstat.utils_1.20-2 base64enc_0.1-3      
[37] pkgconfig_2.0.3       htmltools_0.5.0       parallelly_1.23.0    
[40] fastmap_1.0.1         htmlwidgets_1.5.3     rlang_1.0.2          
[43] shiny_1.5.0           farver_2.0.3          generics_0.1.0       
[46] zoo_1.8-8             jsonlite_1.7.2        ica_1.0-2            
[49] dplyr_1.0.8           magrittr_2.0.1        patchwork_1.1.1      
[52] Matrix_1.3-2          Rcpp_1.0.5            IRkernel_1.1.1       
[55] munsell_0.5.0         fansi_0.4.1           abind_1.4-5          
[58] reticulate_1.18       lifecycle_1.0.1       stringi_1.5.3        
[61] MASS_7.3-53           Rtsne_0.15            plyr_1.8.6           
[64] grid_4.0.3            parallel_4.0.3        listenv_0.8.0        
[67] promises_1.1.1        ggrepel_0.9.0         crayon_1.5.1         
[70] deldir_0.2-3          miniUI_0.1.1.1        lattice_0.20-41      
[73] IRdisplay_0.7.0       cowplot_1.1.1         splines_4.0.3        
[76] tensor_1.5            pillar_1.7.0          igraph_1.3.0         
[79] uuid_0.1-4            future.apply_1.7.0    reshape2_1.4.4       
[82] codetools_0.2-18      leiden_0.3.6          glue_1.6.2           
[85] evaluate_0.14         data.table_1.13.6     png_0.1-7            
[88] vctrs_0.4.0           httpuv_1.5.4          polyclip_1.10-0      
[91] gtable_0.3.0          RANN_2.6.1            purrr_0.3.4          
[94] tidyr_1.2.0           scattermore_0.7       future_1.21.0        
[97] assertthat_0.2.1      ggplot2_3.3.3         mime_0.9             
[100] xtable_1.8-4          RSpectra_0.16-0       later_1.1.0.1        
[103] viridisLite_0.3.0     survival_3.2-7        tibble_3.1.6         
[106] cluster_2.1.0         globals_0.14.0        fitdistrplus_1.1-3   
[109] ellipsis_0.3.2        ROCR_1.0-11          




















