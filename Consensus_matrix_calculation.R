library(doParallel)

consensus_mat <- function(i) { 
  
  setwd("/home/uqdmizik/90days/gwas_test/pipeline")
  meta <- read.table("RC_all_cluster_assignments.txt")
  meta <- as.matrix(meta)
  all <- round(seq(0,nrow(meta), length.out = 30), 0)
  to_do <- seq(all[i]+1, all[i+1], 1)
  m <- matrix(nrow = length(to_do), ncol = nrow(meta))
  genes <- rownames(meta)
  rownames(m) <- genes[to_do]
  colnames(m) <- genes
  
  
  for (x in 1:length(to_do)) {
    for (f in 1:nrow(meta)) {
      
      
      m[x,f] <- length(intersect(meta[rownames(meta) == genes[to_do[x]],], meta[rownames(meta) == genes[f],]))
    }
    print(x) 
  }
  
  
  write.table(m, paste("set", i, "RC_consensus_matrix", sep = "_"))
  
} 


cl <- makeCluster(30)
registerDoParallel(cl)
clusterApply(cl, 1:29, consensus_mat)
stopCluster(cl)