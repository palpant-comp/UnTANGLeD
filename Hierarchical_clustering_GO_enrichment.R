
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cluster)

for (i in 1:29) {
    
    if (i == 1) {
      
      cm <- read.table(paste("set", i, "RC_consensus_matrix", sep = "_"))
      
    } else {
      
      cm <- rbind(cm, read.table(paste("set", i, "RC_consensus_matrix", sep = "_")))
    }
    
    print(i)
    
  }
  
  write.table(cm, paste("RC", "consensus_matrix_full.txt", sep = "_"))
  saveRDS(cm, paste("RC", "consensus_matrix_full.rds", sep = "_"))

  m <- cm
  #### m is matrix ####
  m <- as.matrix(m)
  m <- abs(m-100)
  
  #### m2 is distance matrix ####
  m2 <- as.dist(m)
  
  
  ### clustering with best method ###
  hcl <- hclust(m2, method = "ward.D2")

  sub_grp <- cutree(hcl, k = 242)
  
  #write.table(sub_grp, paste("RC", "cluster_assignments_consensus.txt", sep = "_"))
  
  print(table(sub_grp))
  
  ### individual silhouette measure ###
  sil_d <- silhouette(sub_grp, dist = m2)
  sil_d <- as.data.frame(sil_d[order(sil_d[,1]),])
  
  write.table(sil_d, paste("RC", "silhouette_scores_individual_clus.txt", sep = "_"))
  
  levels <- unique(sub_grp)
  background <- as.character(names(sub_grp))
  length(background)
  GO_df <- as.data.frame(matrix(ncol = 11, nrow = 10000))
  rownames(GO_df) <- 1:10000
  
  ####LOOP TO RUN THROUGH ALL CLUSTERS
  
  
  
  #geneset
  GO_df <- NA
  for (i in 1:length(levels)) {
    cluster <- names(sub_grp[sub_grp == levels[i]])
    head(cluster)
    print(length(cluster))
    
    ### actual GO enrichment
    
    
    ego <- enrichGO(           gene          = cluster,
                               universe      = background,
                               keyType       = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "ALL",
                               pAdjustMethod = "fdr",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
    
    
    
    ### 'simplified version"
    ego <- as.data.frame(ego)
    print(head(ego$Description))

    if (nrow(ego) > 0) {

    ego$Cluster <- levels[i]

    if (is.na(GO_df)) {
      GO_df <- ego
    } else {
      GO_df <- rbind(GO_df, ego)
    } }
    
    
    
  }
  colnames(GO_df) <- c("ONTOLOGY", "ID","Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Cluster")
  GO_df <- na.omit(GO_df)
  write.table(GO_df, "RC_GO_enrichment.txt")
  

