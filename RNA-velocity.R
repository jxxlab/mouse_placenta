
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
pbmc <- subset(pbmc, subset = percent.mt > 1);pbmc
DimPlot(pbmc, label=T,group.by = 'cluster_21',cols=C21)


ldat  <- read.loom.matrices("mpl10.loom")
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures =2000)
x=VariableFeatures(pbmc);length(x)
emat <- ldat$spliced ; nmat <- ldat$unspliced
emat <- emat[x,colnames(pbmc$RNA@counts)]; nmat <- nmat[x,colnames(pbmc$RNA@counts)]
dim(emat) 

cluster.label <- pbmc$cluster_21 ;  cell.colors<- pagoda2:::fac2col(cluster.label,level.colors =C21)  ; 
emb <- pbmc$umap@cell.embeddings ;  dim(emb)
cell.dist <- as.dist(1-armaCor(t(pbmc@reductions$pca@cell.embeddings)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.1)
length(intersect(rownames(emat),rownames(nmat)))  ;    fit.quantile <- 0.02

gene1=intersect(rownames(emat),rownames(nmat))
genex=readLines(con = "murk_genes_C17.csv")
emat=emat[setdiff(gene1,genex),]
nmat=nmat[setdiff(gene1,genex),]
length(intersect(rownames(emat),rownames(nmat)))
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,n.cores=40,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)
show.velocity.on.embedding.cor(emb,rvel.cd,n=900,scale='sqrt',cell.colors=ac(cell.colors,alpha=1),cex=0.6,n.cores=30,  
                                  arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=10,grid.n=30,arrow.lwd=1.5, do.par=T,cell.border.alpha =0)









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
[1] velocyto.R_0.6     pagoda2_0.1.1      igraph_1.2.4.2     Matrix_1.2-18     
[5] RColorBrewer_1.1-2 Seurat_3.1.5      

loaded via a namespace (and not attached):
  [1] tsne_0.1-3          nlme_3.1-139        bitops_1.0-6       
[4] bit64_0.9-7         RcppAnnoy_0.0.14    httr_1.4.2         
[7] repr_1.1.4.9000     sctransform_0.2.1   tools_3.6.0        
[10] utf8_1.1.4          R6_2.4.1            irlba_2.3.3        
[13] KernSmooth_2.23-15  BiocGenerics_0.32.0 mgcv_1.8-28        
[16] uwot_0.1.8          DBI_1.1.0           lazyeval_0.2.2     
[19] colorspace_1.4-1    npsurv_0.4-0        gridExtra_2.3      
[22] tidyselect_1.1.2    bit_4.0.4           compiler_3.6.0     
[25] Biobase_2.46.0      cli_3.2.0           hdf5r_1.3.1        
[28] plotly_4.9.2        labeling_0.3        triebeard_0.3.0    
[31] caTools_1.18.0      scales_1.1.0        lmtest_0.9-37      
[34] ggridges_0.5.2      pbapply_1.4-2       rappdirs_0.3.1     
[37] pbdZMQ_0.3-3        stringr_1.4.0       digest_0.6.24      
[40] base64enc_0.1-3     pkgconfig_2.0.3     htmltools_0.4.0    
[43] htmlwidgets_1.5.1   rlang_1.0.2         farver_2.0.3       
[46] generics_0.0.2      zoo_1.8-7           jsonlite_1.8.0     
[49] ica_1.0-2           gtools_3.8.1        dplyr_1.0.8        
[52] magrittr_2.0.3      patchwork_1.0.1     Rcpp_1.0.3         
[55] IRkernel_1.1        munsell_0.5.0       fansi_0.4.1        
[58] ape_5.3             reticulate_1.14     lifecycle_1.0.1    
[61] stringi_1.4.5       MASS_7.3-51.4       gplots_3.0.1.2     
[64] Rtsne_0.15          plyr_1.8.5          grid_3.6.0         
[67] parallel_3.6.0      gdata_2.18.0        listenv_0.8.0      
[70] ggrepel_0.8.1       crayon_1.5.1        lattice_0.20-38    
[73] IRdisplay_0.7.0     cowplot_1.0.0       splines_3.6.0      
[76] pillar_1.7.0        uuid_0.1-4          rjson_0.2.20       
[79] future.apply_1.4.0  reshape2_1.4.3      codetools_0.2-16   
[82] leiden_0.3.5        glue_1.6.2          evaluate_0.14      
[85] lsei_1.2-0          pcaMethods_1.78.0   data.table_1.14.2  
[88] urltools_1.7.3      vctrs_0.4.0         png_0.1-7          
[91] gtable_0.3.0        RANN_2.6.1          purrr_0.3.4        
[94] tidyr_1.2.0         future_1.16.0       assertthat_0.2.1   
[97] ggplot2_3.3.3       rsvd_1.0.2          survival_3.1-8     
[100] viridisLite_0.3.0   tibble_3.1.6        Rook_1.1-1         
[103] cluster_2.0.8       globals_0.12.5      brew_1.0-6         
[106] fitdistrplus_1.0-14 ellipsis_0.3.2      ROCR_1.0-7    
