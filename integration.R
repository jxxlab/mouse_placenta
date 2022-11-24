library(Seurat)
library(RColorBrewer)
library(SeuratDisk)

pbmc <- readRDS("/home/jxx/MPL/nu16836.rds")
pbmc$cluster_17E=Idents(pbmc)
Idents(pbmc)="cluster_17E"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI")) <- "21"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII")) <- "20"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII Precursor")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT Precursor")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "Glycogen Cells")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 2")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 1")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC Precursor")) <- "6-7-8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI Precursor")) <- "4-5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP 2")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP")) <- "2"
pbmc$cluster_21=Idents(pbmc)

Idents(pbmc)="cluster_17E"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI")) <- "21 SynTI"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII")) <- "20 SynTII"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII Precursor")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT Precursor")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "Glycogen Cells")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 2")) <- "14-15   SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 1")) <- "14-15   SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC")) <- "9   S-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC Precursor")) <- "6-7-8 S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI Precursor")) <- "4-5 SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP 2")) <- "3   LaTP 2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP")) <- "2   LaTP"
pbmc$type_21=Idents(pbmc)

Idents(pbmc)='GA'
table(Idents(pbmc))
new.cluster.ids <- c("E09.5n","E10.5n",'E12.5n','E14.5n')
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
table(Idents(pbmc))
pbmc$stage=Idents(pbmc)
pbmc$sample=Idents(pbmc)
nuclei <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data, min.cells = 3);nuclei

pbmc <- readRDS("MTR15682.rds");pbmc
cells <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data, min.cells = 3);cells

pbmc <- LoadH5Seurat("/sdb/jxx/MTR/V3MTR8067.h5Seurat");pbmc
Idents(pbmc)="stage"
new.cluster.ids <- c('E07.5v3','E08.5v3',"E09.5v3",'E11.5v3', 'E13.5v3')
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
pbmc$stage=Idents(pbmc)
pbmc$sample=Idents(pbmc)
v3 <- CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data, min.cells = 3);v3



nuclei$origin="Nuclei"
cells$origin="Cells"
v3$origin="v3"
pbmc <- merge(x =cells, y =c(nuclei,v3) , project = "MTR")
pbmc;table(pbmc$origin)
strn=rownames(nuclei)
strc=rownames(cells)
strv3=rownames(v3)
str=intersect(intersect(strn,strc),strv3);length(str);head(str)
count=(pbmc$RNA@counts[str,]);dim(count)
pbmc <- CreateSeuratObject(count, meta.data = pbmc@meta.data, min.cells = 3, project = "MTR");pbmc

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(pbmc))
pancreas=pbmc
pancreas.list <- SplitObject(object = pancreas, split.by = "stage")
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = 8000, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = 30, verbose = FALSE)
DimPlot(pancreas.integrated);  ElbowPlot(pancreas.integrated)
pbmc <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:12,verbose=F)
DimPlot(pbmc, label=T)

Idents(pbmc)="cluster_17E"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P3")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P2")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1")) <- "P1-2-3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT Precursor")) <- "17"
Idents(pbmc, cells = WhichCells(pbmc,idents = "Glycogen Cells")) <- "16"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 2")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 1")) <- "14-15"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC")) <- "9"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC Precursor")) <- "6-7-8"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI")) <- "21"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII")) <- "20"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI Precursor")) <- "4-5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP 2")) <- "3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII Precursor")) <- "19"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP")) <- "2"
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
Idents(pbmc)="cluster_17E"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P3")) <- "P1-2-3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P2")) <- "P1-2-3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "P1")) <- "P1-2-3 Biopotential progenitor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E2")) <- "E2 The progenitor of S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E1")) <- "E1 The progenitor of SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SpT Precursor")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "Glycogen Cells")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 2")) <- "14-15   SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "JZP 1")) <- "14-15   SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC")) <- "9   S-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "S-TGC Precursor")) <- "6-7-8 S-TGC Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI")) <- "21 SynTI"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII")) <- "20 SynTII"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTI Precursor")) <- "4-5 SynTI Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP 2")) <- "3   LaTP 2"
Idents(pbmc, cells = WhichCells(pbmc,idents = "SynTII Precursor")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "LaTP")) <- "2   LaTP"
Idents(pbmc, cells = WhichCells(pbmc,idents = "19")) <- "19 SynTII Precursor"
Idents(pbmc, cells = WhichCells(pbmc,idents = "18")) <- "18 EPC Migratory Cell"
Idents(pbmc, cells = WhichCells(pbmc,idents = "17")) <- "17 SpA-TGC"
Idents(pbmc, cells = WhichCells(pbmc,idents = "16")) <- "16 Gly-T"
Idents(pbmc, cells = WhichCells(pbmc,idents = "15")) <- "14-15   SpT"
Idents(pbmc, cells = WhichCells(pbmc,idents = "14")) <- "14-15   SpT"
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

C21=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1:4,7,9,22:24,12:14,27,30,29,5,15,10,18)]
Idents(pbmc)="origin"
new.cluster.ids <- c('Cells_V2','Nuclei_V3','Cells_V3')
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids, reorder.numeric = T)
Idents(pbmc, cells = WhichCells(pbmc,idents = "Cells_V3")) <- "Cells_V3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "Cells_V2")) <- "Cells_V2"
pbmc$origin=Idents(pbmc)
table(Idents(pbmc))
Idents(pbmc)='cluster_21'
DimPlot(pbmc, label=T,split.by="origin",cols=C21)
Idents(pbmc)='sample'
Idents(pbmc, cells = WhichCells(pbmc,idents = "E13.5v3")) <- "E13.5v3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E13.5")) <- "E13.5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E12.5n")) <- "E12.5n"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E12.5h")) <- "E12.5h"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E12.5")) <- "E12.5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E11.5v3")) <- "E11.5v3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E11.5h")) <- "E11.5h"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E10.5n")) <- "E10.5n"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E10.5h")) <- "E10.5h"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E10.5")) <- "E10.5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E09.5n")) <- "E9.5n"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E09.5v3")) <- "E9.5v3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E09.5")) <- "E9.5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E08.5v3")) <- "E8.5v3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E08.5")) <- "E8.5"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E07.5v3")) <- "E7.5v3"
Idents(pbmc, cells = WhichCells(pbmc,idents = "E07.5")) <- "E7.5"
DimPlot(pbmc, label=F,cols=colorRampPalette(rainbow(9)[c(4,3,6,7,2,1,9,8)])(19))




session_info()

R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /usr/lib64/R/lib/libRblas.so
LAPACK: /usr/lib64/R/lib/libRlapack.so

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
[1] SeuratDisk_0.0.0.9013 RColorBrewer_1.1-2    Seurat_3.1.5         

loaded via a namespace (and not attached):
[1] tsne_0.1-3          nlme_3.1-139        bitops_1.0-6       
[4] bit64_0.9-7         RcppAnnoy_0.0.14    httr_1.4.2         
[7] repr_1.1.4.9000     sctransform_0.2.1   tools_3.6.0        
[10] utf8_1.1.4          R6_2.4.1            irlba_2.3.3        
[13] KernSmooth_2.23-15  uwot_0.1.8          DBI_1.1.0          
[16] lazyeval_0.2.2      colorspace_1.4-1    withr_2.5.0        
[19] npsurv_0.4-0        gridExtra_2.3       tidyselect_1.1.2   
[22] bit_4.0.4           compiler_3.6.0      cli_3.2.0          
[25] hdf5r_1.3.1         plotly_4.9.2        labeling_0.3       
[28] caTools_1.18.0      scales_1.1.0        lmtest_0.9-37      
[31] ggridges_0.5.2      pbapply_1.4-2       pbdZMQ_0.3-3       
[34] stringr_1.4.0       digest_0.6.24       base64enc_0.1-3    
[37] pkgconfig_2.0.3     htmltools_0.4.0     htmlwidgets_1.5.1  
[40] rlang_1.0.2         farver_2.0.3        generics_0.0.2     
[43] zoo_1.8-7           jsonlite_1.8.0      ica_1.0-2          
[46] gtools_3.8.1        dplyr_1.0.8         magrittr_2.0.3     
[49] patchwork_1.0.1     Matrix_1.2-18       Rcpp_1.0.3         
[52] IRkernel_1.1        munsell_0.5.0       fansi_0.4.1        
[55] ape_5.3             reticulate_1.14     lifecycle_1.0.1    
[58] stringi_1.4.5       MASS_7.3-51.4       gplots_3.0.1.2     
[61] Rtsne_0.15          plyr_1.8.5          grid_3.6.0         
[64] parallel_3.6.0      gdata_2.18.0        listenv_0.8.0      
[67] ggrepel_0.8.1       crayon_1.5.1        lattice_0.20-38    
[70] IRdisplay_0.7.0     cowplot_1.0.0       splines_3.6.0      
[73] pillar_1.7.0        igraph_1.2.4.2      uuid_0.1-4         
[76] future.apply_1.4.0  reshape2_1.4.3      codetools_0.2-16   
[79] leiden_0.3.5        glue_1.6.2          evaluate_0.14      
[82] lsei_1.2-0          data.table_1.14.2   vctrs_0.4.0        
[85] png_0.1-7           gtable_0.3.0        RANN_2.6.1         
[88] purrr_0.3.4         tidyr_1.2.0         future_1.16.0      
[91] assertthat_0.2.1    ggplot2_3.3.3       rsvd_1.0.2         
[94] RSpectra_0.16-0     survival_3.1-8      viridisLite_0.3.0  
[97] tibble_3.1.6        cluster_2.0.8       globals_0.12.5     
[100] fitdistrplus_1.0-14 ellipsis_0.3.2      ROCR_1.0-7 



