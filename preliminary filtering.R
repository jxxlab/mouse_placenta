
library(Seurat)
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices_mex/GRCm38.89")
pbmc <- CreateSeuratObject(counts = pbmc.data, min.features = 200, names.field = 2, names.delim = "-",  project = "MTR")
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")
pbmc <- subset(x = pbmc, subset =nFeature_RNA > 500 & percent.mt<10);pbmc
table(Idents(pbmc))

pbmc2=pbmc

pbmc <- SubsetData(object = pbmc2, ident.use = 1)
write.csv(pbmc@meta.data, "matrix/dE07.5.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 2)
write.csv(pbmc@meta.data, "matrix/dE08.5.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 3)
write.csv(pbmc@meta.data, "matrix/dE09.5.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 5)
write.csv(pbmc@meta.data, "matrix/dE10.5.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 6)
write.csv(pbmc@meta.data, "matrix/dE10.5h.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 7)
write.csv(pbmc@meta.data, "matrix/dE11.5h.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 8)
write.csv(pbmc@meta.data, "matrix/dE12.5h.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 9)
write.csv(pbmc@meta.data, "matrix/dE12.5.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 10)
write.csv(pbmc@meta.data, "matrix/dE13.5.csv")
pbmc <- SubsetData(object = pbmc2, ident.use = 11)
write.csv(pbmc@meta.data, "matrix/dE14.5.csv")

pbmc=pbmc2
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =1)])
write.csv(pbmc.data, "matrix/E07.5.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =2)])
write.csv(pbmc.data, "matrix/E08.5.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =3)])
write.csv(pbmc.data, "matrix/E09.5.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =5)])
write.csv(pbmc.data, "matrix/E10.5.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =6)])
write.csv(pbmc.data, "matrix/E10.5h.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =7)])
write.csv(pbmc.data, "matrix/E11.5h.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =8)])
write.csv(pbmc.data, "matrix/E12.5.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =9)])
write.csv(pbmc.data, "matrix/E12.5h.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =10)])
write.csv(pbmc.data, "matrix/E13.5.csv")
pbmc.data <- as.matrix(x = pbmc$RNA@counts[, WhichCells(pbmc, ident =11)])
write.csv(pbmc.data, "matrix/E14.5.csv")


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
[1] Seurat_3.1.5

loaded via a namespace (and not attached):
[1] tsne_0.1-3          nlme_3.1-139        bitops_1.0-6        fs_1.3.1            usethis_1.5.0       devtools_2.0.2     
[7] RcppAnnoy_0.0.12    RColorBrewer_1.1-2  httr_1.4.0          rprojroot_2.0.2     sctransform_0.2.0   tools_3.6.0        
[13] R6_2.4.0            irlba_2.3.3         KernSmooth_2.23-15  uwot_0.1.8          lazyeval_0.2.2      colorspace_1.4-1   
[19] npsurv_0.4-0        withr_2.1.2         gridExtra_2.3       tidyselect_1.1.0    prettyunits_1.0.2   processx_3.3.1     
[25] compiler_3.6.0      cli_1.1.0           desc_1.2.0          plotly_4.9.0        caTools_1.17.1.2    scales_1.0.0       
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
