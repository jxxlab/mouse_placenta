library(monocle)
library(Seurat)
library(RColorBrewer)
library(pagoda2)
library(velocyto.R)

pbmc=readRDS("MTR15682.rds");pbmc
write.table(colnames(pbmc), file="cell15682.txt",row.names = F, col.names = F,quote = F)  
# modify the cell name with excel , and make the cell name be the same with the cell name produced by velocyto
cell=readLines("cell15682x.txt");length(cell)
matrix=pbmc$RNA@counts
colnames(matrix)=cell
meta=pbmc@meta.data
rownames(meta)=cell

pbmc <- CreateSeuratObject(matrix, meta.data =meta, min.cells = 3, project = "mtr");pbmc
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(pbmc))
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc,verbose = F)
DimPlot(pbmc, reduction = "pca", label=T)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:11,umap.method ='umap-learn',metric= 'correlation')
C21=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1:4,7,9,22:24,12:14,27,30,15,10,18)]
DimPlot(pbmc, label=T,group.by = 'cluster_21',cols=C21)

umap=readRDS("MTR15682.rds")$umap@cell.embeddings
rownames(umap)=cell
pbmc$umap@cell.embeddings=umap
DimPlot(pbmc, label=T,group.by = 'cluster_21',cols=C21)

Idents(pbmc)="cluster_17E"
pbmc <- subset(pbmc, idents =c(1,"P1","P2","P3","E2",6:9));pbmc
C9=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1,6:9,10,17,18,20)]
C5=colorRampPalette(brewer.pal(12, "Paired"))(30)[c(1,7,9,10,18)]
DimPlot(pbmc,label=T,label.size = 8,group.by="cluster_17E" ,cols=C9 )
DimPlot(pbmc,label=T,label.size = 8,group.by="cluster_21",cols=C5  )

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 5), dispersion.cutoff = c(0.5, Inf))
length(VariableFeatures(pbmc))
pbmc <- ScaleData(pbmc); pbmc <- RunPCA(pbmc,verbose = F)
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc) 
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:6,umap.method ='umap-learn',metric= 'correlation')
options(repr.plot.width=7.2, repr.plot.height=6)
DimPlot(pbmc,label=T,label.size = 8,group.by="cluster_17E" ,cols=C9  )
DimPlot(pbmc,label=T,label.size = 8,group.by="cluster_21",cols=C5  )


newimportCDS <- function (otherCDS, seurat_scale=F, import_all = FALSE) 
{
  if (class(otherCDS)[1] == "Seurat") {
    requireNamespace("Seurat")
    if (!seurat_scale) {
      data <- otherCDS@assays$RNA@counts
    } else {
      data <- otherCDS@assays$RNA@scale.data
    }
    if (class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data), 
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    #lowerDetectionLimit <- otherCDS@is.expr
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
      expr <- "negbinomial.size"
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
      expr <- "unimormal"
    }
    else {
      expressionFamily <- tobit()
      expr <- "tobit"
    }
    print(paste0("expressionFamily ",expr))
    # valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  #lowerDetectionLimit = lowerDetectionLimit,
                                  expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }
      else {
        mist_list <- otherCDS
      }
    }
    else {
      mist_list <- list()
    }
    if ("var.genes" %in% slotNames(otherCDS)) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
  }
  else if (class(otherCDS)[1] == "SCESet") {
    requireNamespace("scater")
    message("Converting the exprs data in log scale back to original scale ...")
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if ("is.expr" %in% slotNames(otherCDS)) 
      lowerDetectionLimit <- otherCDS@is.expr
    else lowerDetectionLimit <- 1
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }
    else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }
    else {
      expressionFamily <- tobit()
    }
    if (import_all) {
      mist_list <- otherCDS
    }
    else {
      mist_list <- list()
    }
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                                  lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    monocle_cds@auxOrderingData$scran <- mist_list
  }
  else {
    stop("the object type you want to export to is not supported yet")
  }
  return(monocle_cds)
}


HSMM <- newimportCDS(CreateSeuratObject(pbmc$RNA@counts, meta.data = pbmc@meta.data, min.cells = 3, project = "mtr")  , import_all = TRUE);HSMM
HSMM <- estimateSizeFactors(HSMM);HSMM <- estimateDispersions(HSMM)
HSMM_myo <-HSMM
disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table, mean_expression >= 0.01& dispersion_empirical >= 1 * dispersion_fit)$gene_id
length(ordering_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM_myo <-  reduceDimension(HSMM_myo, method = 'DDRTree')
HSMM_myo <-  orderCells(HSMM_myo, reverse =F)
plot_cell_trajectory(HSMM_myo, color_by = "cluster_17E")+scale_colour_manual(values = C9)
plot_cell_trajectory(HSMM_myo, color_by = "cluster_21")+scale_colour_manual(values = C5)
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
plot_cell_trajectory(HSMM_myo, color_by ="cluster_21")+scale_colour_manual(values = C5) +facet_wrap(~stage, ncol = 4)
plot_genes_in_pseudotime(HSMM_myo[c("Hand1","Bhlhe40","Cdx1", "Mta3")], relative_expr = F, ncol = 4,trend_formula = "~ sm.ns(Pseudotime, df=5)", color_by = "stage") + 
  scale_colour_manual(values = rainbow(9)[c(4,3,6,7,2,1,9,8)])

ldat  <- read.loom.matrices("mpl10.loom")
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures =2000)
x=VariableFeatures(pbmc);length(x)
emat <- ldat$spliced ; nmat <- ldat$unspliced
emat <- emat[x,colnames(pbmc$RNA@counts)]; nmat <- nmat[x,colnames(pbmc$RNA@counts)]
dim(emat) 
cluster.label <- pbmc$cluster_17E ;  cell.colors<- pagoda2:::fac2col(cluster.label,level.colors =C9)  ; 
emb <- pbmc$umap@cell.embeddings ;  dim(emb)
cell.dist <- as.dist(1-armaCor(t(pbmc@reductions$pca@cell.embeddings)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.1)
length(intersect(rownames(emat),rownames(nmat)))  ;    fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,n.cores=30,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',  cell.colors=ac(cell.colors,alpha=1),cex=0.4,n.cores=30,
                               arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=15,arrow.lwd=1.5, do.par=F,cell.border.alpha = 0.1)

emb <- t(HSMM_myo@reducedDimS);  dim(emb)
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',  cell.colors=ac(cell.colors,alpha=1),cex=0.4,
                               arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=15,arrow.lwd=1.5, do.par=F,cell.border.alpha = 0.1)



sessionInfo()


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
[1] splines   stats4    parallel  stats     graphics  grDevices utils    
[8] datasets  methods   base     

other attached packages:
[1] velocyto.R_0.6      pagoda2_0.1.1       igraph_1.2.4.2     
[4] RColorBrewer_1.1-2  Seurat_3.1.5        monocle_2.14.0     
[7] DDRTree_0.1.5       irlba_2.3.3         VGAM_1.1-2         
[10] ggplot2_3.3.3       Biobase_2.46.0      BiocGenerics_0.32.0
[13] Matrix_1.2-18      

loaded via a namespace (and not attached):
  [1] Rtsne_0.15           colorspace_1.4-1     rjson_0.2.20        
[4] ellipsis_0.3.2       ggridges_0.5.2       IRdisplay_0.7.0     
[7] base64enc_0.1-3      proxy_0.4-23         farver_2.0.3        
[10] leiden_0.3.5         listenv_0.8.0        npsurv_0.4-0        
[13] urltools_1.7.3       bit64_0.9-7          ggrepel_0.8.1       
[16] fansi_0.4.1          codetools_0.2-16     docopt_0.6.1        
[19] lsei_1.2-0           IRkernel_1.1         jsonlite_1.8.0      
[22] ica_1.0-2            cluster_2.0.8        png_0.1-7           
[25] pheatmap_1.0.12      uwot_0.1.8           sctransform_0.2.1   
[28] compiler_3.6.0       httr_1.4.2           assertthat_0.2.1    
[31] lazyeval_0.2.2       limma_3.42.2         cli_3.2.0           
[34] htmltools_0.4.0      tools_3.6.0          rsvd_1.0.2          
[37] gtable_0.3.0         glue_1.6.2           RANN_2.6.1          
[40] reshape2_1.4.3       dplyr_1.0.8          rappdirs_0.3.1      
[43] Rcpp_1.0.3           slam_0.1-47          vctrs_0.4.0         
[46] gdata_2.18.0         ape_5.3              nlme_3.1-139        
[49] lmtest_0.9-37        stringr_1.4.0        globals_0.12.5      
[52] lifecycle_1.0.1      gtools_3.8.1         future_1.16.0       
[55] MASS_7.3-51.4        zoo_1.8-7            scales_1.1.0        
[58] pcaMethods_1.78.0    reticulate_1.14      pbapply_1.4-2       
[61] gridExtra_2.3        triebeard_0.3.0      fastICA_1.2-2       
[64] stringi_1.4.5        Rook_1.1-1           caTools_1.18.0      
[67] densityClust_0.3     repr_1.1.4.9000      rlang_1.0.2         
[70] pkgconfig_2.0.3      matrixStats_0.55.0   bitops_1.0-6        
[73] qlcMatrix_0.9.7      evaluate_0.14        lattice_0.20-38     
[76] ROCR_1.0-7           purrr_0.3.4          labeling_0.3        
[79] patchwork_1.0.1      htmlwidgets_1.5.1    bit_4.0.4           
[82] cowplot_1.0.0        tidyselect_1.1.2     RcppAnnoy_0.0.14    
[85] plyr_1.8.5           magrittr_2.0.3       R6_2.4.1            
[88] gplots_3.0.1.2       generics_0.0.2       combinat_0.0-8      
[91] pbdZMQ_0.3-3         DBI_1.1.0            mgcv_1.8-28         
[94] pillar_1.7.0         withr_2.5.0          fitdistrplus_1.0-14 
[97] survival_3.1-8       tibble_3.1.6         future.apply_1.4.0  
[100] tsne_0.1-3           hdf5r_1.3.1          crayon_1.5.1        
[103] uuid_0.1-4           KernSmooth_2.23-15   utf8_1.1.4          
[106] plotly_4.9.2         viridis_0.5.1        grid_3.6.0          
[109] data.table_1.14.2    FNN_1.1.3            HSMMSingleCell_1.6.0
[112] sparsesvd_0.2        digest_0.6.24        tidyr_1.2.0         
[115] brew_1.0-6           munsell_0.5.0        viridisLite_0.3.0   





