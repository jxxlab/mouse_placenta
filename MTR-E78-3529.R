library(Seurat);library("RColorBrewer")

pbmc <- readRDS("MTR15682.rds");pbmc
Idents(pbmc) <- "stage"
pbmc <- subset(pbmc, idents =  c("E07.5","E08.5"))

pbmc <- CreateSeuratObject(pbmc$RNA@counts, meta.data =pbmc@meta.data, min.cells = 3, project = "mtr");pbmc
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(pbmc))
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,verbose = F)
DimPlot(pbmc, reduction = "pca", label=T)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:7,umap.method ='umap-learn',metric= 'correlation')
DimPlot(pbmc, label=T) 


pbmc <- FindNeighbors(object = pbmc, dims = 1:7)
pbmc <- FindClusters(object = pbmc, algorithm = 1, resolution =0.25)
DimPlot(pbmc,  label=T,label.size = 8)


C18=colnames(subset(pbmc,idents=6))
write.table(C18, file="cluster_18_cell.txt",row.names = F, col.names = F,quote = F) 

Idents(pbmc) <- "cluster_19"
Idents(pbmc, cells = C18) <- "18"
Idents(pbmc, cells = WhichCells(pbmc,idents = c(6,14))) <- "13"
Idents(pbmc, cells = WhichCells(pbmc,idents = "12")) <- "12"
Idents(pbmc, cells = WhichCells(pbmc,idents = "11")) <- "11"
Idents(pbmc, cells = WhichCells(pbmc,idents = "10")) <- "10"
Idents(pbmc, cells = WhichCells(pbmc,idents = "5")) <- "4"
Idents(pbmc, cells = WhichCells(pbmc,idents = c("2","19"))) <- "1"
pbmc$cluster_7=Idents(pbmc)
options(repr.plot.width=6.9, repr.plot.height=6)
C7=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1,4,22:24,18,27)]
DimPlot(pbmc, label=T,label.size = 8,cols=C7)

saveRDS(pbmc, file = "MTR-E78-3529.rds")




sessionInfo()

# the output of sessionInfo()
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
[1] RColorBrewer_1.1-2 Seurat_3.1.5      

loaded via a namespace (and not attached):
[1] tsne_0.1-3          nlme_3.1-139        bitops_1.0-6       
[4] RcppAnnoy_0.0.14    httr_1.4.2          repr_1.1.4.9000    
[7] sctransform_0.2.1   tools_3.6.0         utf8_1.1.4         
[10] R6_2.4.1            irlba_2.3.3         KernSmooth_2.23-15 
[13] uwot_0.1.8          DBI_1.1.0           lazyeval_0.2.2     
[16] colorspace_1.4-1    npsurv_0.4-0        gridExtra_2.3      
[19] tidyselect_1.1.2    compiler_3.6.0      cli_3.2.0          
[22] plotly_4.9.2        labeling_0.3        caTools_1.18.0     
[25] scales_1.1.0        lmtest_0.9-37       ggridges_0.5.2     
[28] pbapply_1.4-2       rappdirs_0.3.1      pbdZMQ_0.3-3       
[31] stringr_1.4.0       digest_0.6.24       base64enc_0.1-3    
[34] pkgconfig_2.0.3     htmltools_0.4.0     htmlwidgets_1.5.1  
[37] rlang_1.0.2         farver_2.0.3        generics_0.0.2     
[40] zoo_1.8-7           jsonlite_1.8.0      ica_1.0-2          
[43] gtools_3.8.1        dplyr_1.0.8         magrittr_2.0.3     
[46] patchwork_1.0.1     Matrix_1.2-18       Rcpp_1.0.3         
[49] IRkernel_1.1        munsell_0.5.0       fansi_0.4.1        
[52] ape_5.3             reticulate_1.14     lifecycle_1.0.1    
[55] stringi_1.4.5       MASS_7.3-51.4       gplots_3.0.1.2     
[58] Rtsne_0.15          plyr_1.8.5          grid_3.6.0         
[61] parallel_3.6.0      gdata_2.18.0        listenv_0.8.0      
[64] ggrepel_0.8.1       crayon_1.5.1        lattice_0.20-38    
[67] IRdisplay_0.7.0     cowplot_1.0.0       splines_3.6.0      
[70] pillar_1.7.0        igraph_1.2.4.2      uuid_0.1-4         
[73] future.apply_1.4.0  reshape2_1.4.3      codetools_0.2-16   
[76] leiden_0.3.5        glue_1.6.2          evaluate_0.14      
[79] lsei_1.2-0          data.table_1.14.2   vctrs_0.4.0        
[82] png_0.1-7           gtable_0.3.0        RANN_2.6.1         
[85] purrr_0.3.4         tidyr_1.2.0         future_1.16.0      
[88] assertthat_0.2.1    ggplot2_3.3.3       rsvd_1.0.2         
[91] survival_3.1-8      viridisLite_0.3.0   tibble_3.1.6       
[94] cluster_2.0.8       globals_0.12.5      fitdistrplus_1.0-14
[97] ellipsis_0.3.2      ROCR_1.0-7 









