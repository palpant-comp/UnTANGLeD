#### load packages ####
  library(Seurat)
  library(stringr)
  library(ggplot2)
  library(data.table)
  library(tibble)
  setwd("/home/uqdmizik/90days/gwas_test/pipeline")

  full_data <- "/home/uqdmizik/90days/gwas_test/pipeline/"
  multixcan <- "/home/uqdmizik/90days/gwas_test/multixcan/REAL_PHEN_OUTPUTS/multixcan/"

  tog <- readRDS(paste0(full_data, "full_matrix.rds"))


  t <- tog
  
  t3 <- as.matrix(t)
  t3[is.na(t3)] <- 1
  
  ### remove genes with fewer than one significant association ###
  ### create function to check # of significant associations
  significance <- function(x) {table(x < 10^-4)}
  
  ##### apply significance function to genes
  t <- (apply(t3, 2, significance))
  
  vt<-NULL
  vf<-NULL
  for (i in 1:length(t)){
    vt[i]<-t[[i]][1]
    vf[i]<-t[[i]][2]
  }
  gene.significance <- cbind(vt,vf)
  rownames(gene.significance) <- names(t)
  gene.significance <- gene.significance[order(gene.significance[,2], decreasing = T),]
  gene.significance[is.na(gene.significance)] <- 0
  
  sigg <- rownames(subset(gene.significance, gene.significance[,2] > 1))
  
  t3 <- t3[,colnames(t3) %in% sigg]
  
  t3 <- qchisq(t3, 1, lower.tail = F)
  
  ### replace infinity values ###
  range(t3, finite = TRUE)
  table(is.infinite(t3))
  ## set infinite values to value slightly higher than maximum
  t3[is.infinite(t3)] <- (max(t3[is.finite(t3)]) + 5)
  range(t3, finite = TRUE)
  table(is.infinite(t3))
  
  rownames(t3)[is.na(rownames(t3))] <- "NA.x"
  
  t3 <- t(t3)
  ### Normalise phenotypes
  ### Create Seurat Object ###
  Sobject <- CreateSeuratObject(counts = t3, min.cells = 0, min.features  = 0, project = "genometrial")
  
  ### RC normalization ##
  t6 <- NormalizeData(object = Sobject, normalization.method = "RC", margin = 1, scale.factor = 1000)
  g <- t(as.matrix(GetAssayData(t6, slot = 'data')))
  
  # Recreate Seurat object with RC-normalized data (g)
  Sobject <- CreateSeuratObject(counts = g, min.cells = 0, min.features  = 0, project = "genometrial")
  
  # Redundant: this log-normalizes, but it's overwritten immediately after
  # t6 <- NormalizeData(object = Sobject, normalization.method = "LogNormalize", margin = 1, scale.factor = 1000)
  
  t6 <- Sobject  # Continue with object holding raw/RC-normalized matrix
  t6 <- SetAssayData(t6, slot = 'data', g)
  
  ### Find Variable Features (Not used later) ###
  # t6 <- FindVariableFeatures(t6, selection.method = "vst", mean.cutoff = c(0.1, 10), dispersion.cutoff = c(0.5,20), nfeatures = 1394)
  # head(VariableFeatures(t6))
  # length(VariableFeatures(t6))
  
  ## Scale Data (Overwritten immediately) ##
  # t6 <- ScaleData(object = t6)
  t6 <- SetAssayData(object = t6, new.data = g, slot = 'scale.data')
  
  # PCA
  t7 <- RunPCA(object = t6, verbose = F)
  p1 <- DimPlot(object = t7, dims = c(1,2), reduction = "pca") 
  p1

  ## ACTUAL Clustering ##
  t8 <- FindNeighbors(t7, reduction = "pca", dims = 1:10)
  
  res <- seq(0.2, 20, 0.2)

    #### loop resolutions ####
    for(v in 1:100) {
      t8 <- FindClusters(t8, resolution = res[v], n.start = 25)
      
      
      meta <- t8@meta.data
      head(meta)
      
      if (v == 1) {
        meta <- meta[order(rownames(meta)),]
        final_meta <- cbind(rownames(meta), meta$seurat_clusters)        
      } else {
        meta <- meta[order(rownames(meta)),]
        final_meta <- cbind(final_meta, meta$seurat_clusters)
        
      }
    }
  
  rownames(final_meta) <- final_meta[,1]
  final_meta <- final_meta[,2:101]
  colnames(final_meta) <- res
  
  for (v in 1:100) {
    final_meta[,v] <- paste(final_meta[,v], res[v], sep = "_")
  }
  
  write.table(final_meta, "RC_all_cluster_assignments.txt", row.names = T, col.names = T)
  


