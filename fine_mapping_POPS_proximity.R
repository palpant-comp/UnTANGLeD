### fine mapping - POPS method for proximity based genes ###

library(cluster)
library(stringr)
library(ggplot2)



  todo <- list.files(pattern="*clump2.clumped")
  codes <- str_split_fixed(todo, "_", 2)[,1]
  genes <- read.table("glist-hg38", header = F)
  colnames(genes) <- c("CHR", "start", "end", "gene")
  full_data <- "/home/uqdmizik/90days/gwas_test/pipeline/"
  multixcan <- "/home/uqdmizik/90days/gwas_test/multixcan/REAL_PHEN_OUTPUTS/multixcan/"

  tog <- readRDS(paste0(full_data, "full_matrix.rds"))

  phenotypes <- list.files(multixcan)
  newpcodes <- c(1, 4, 4, 8, 8, 10, 12, 12, 15, 18, 18, 20, 6)
  phenotypes <- phenotypes[newpcodes]
  inflection = 242

  	setwd("/home/uqdmizik/90days/gwas_test/pipeline/phenotype_correlation")


  cormat <- readRDS(paste("base", "gene_correlation.rds", sep = "_"))

  	setwd("/home/uqdmizik/90days/gwas_test/pipeline")


  #cm <- readRDS(paste(pcodes[13], "consensus_matrix_full.rds", sep = "_"))
  cm <- readRDS("RC_consensus_matrix.rds")


rownames(cormat) <- colnames(cm)
colnames(cormat) <- colnames(cm)

m <- cm
	#### m is matrix ####
	m <- as.matrix(m)
	m <- abs(m-100)

	#### m2 is distance matrix ####
	m2 <- as.dist(m)


	### clustering with best method ###
	hcl <- hclust(m2, method = "ward.D2")
	sub_grp <- cutree(hcl, k = inflection)

	sil_d <- silhouette(sub_grp, dist = m2)
	sil_d <- as.data.frame(sil_d[order(sil_d[,1]),])

	library("biomaRt")

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


results <- list()
undetected <- list()
phenotypes <- unique(phenotypes)
phenotypes <- c(phenotypes[1:2], phenotypes[9], phenotypes[3:8])
setwd("/scratch/90days/uqdmizik/gwas_test/fine_mapping")


for (x in 1:length(codes)) {

	
	fm <- read.table(paste(codes[x], "loci_genes.txt", sep = "_"))

	full <- fm$genes %>% paste(collapse = "___") %>% str_split("___") %>% unlist() %>% unique()

	gene_names <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
      filters = 'hgnc_symbol', 
      values = full, 
      mart = ensembl,
      useCache = FALSE)

	fmcomp <- data.frame(loci = fm$SNP, pval = fm$P, ngenes = NA, 
		ngenes_e = NA, clusters = NA, genes_pn = NA, genes_p = NA, genes_np = NA )

	chr <- unique(fm$CHR)
		fm$genes_n <- NA
		fm$genes_pn <- NA
		fm$genes_p <- NA
		fm$genes_total <- NA

	for (i in 1:length(chr)) {

		leftout <- fm[fm$CHR == chr[i],]$genes %>% paste(collapse = "___") %>% str_split("___") %>% unlist()
		leftout <- gene_names[gene_names$hgnc_symbol %in% leftout, ]

		remaining <- fm[!(fm$CHR == chr[i]),]$genes %>% str_split("___")

		pval <- vector()
		or <- vector()

		for (f in 1:242) {
			clus <- names(sub_grp[sub_grp == f])

			counter <- 0
			for (z in 1:length(remaining)) {

				temp <- gene_names[gene_names$hgnc_symbol %in% remaining[[z]],]

				if (length(intersect(clus, temp$ensembl_gene_id)) > 0) {
					counter <- counter + 1
				} 

			}




			fish <- fisher.test(matrix(c(counter,
					(length(clus) - counter), 
				   (length(remaining) - counter),
				   length(setdiff(names(sub_grp), clus)) - length(remaining)), nrow = 2, ncol = 2))
			pval[f] <- fish$p.value
			or[f] <- fish$estimate
		}

		names(pval) <- 1:242
		names(or) <- 1:242

		asclus <- names(pval[pval < 0.01])
		asclus <- names(or[(names(or) %in% asclus) & (or > 1)])

		print(asclus)


		fmgenes <- names(sub_grp[sub_grp %in% asclus])

		leftout <- leftout[leftout$ensembl_gene_id %in% names(sub_grp),]

		for (f in 1:nrow(fm[fm$CHR == chr[i],])) {

			print(f)

			genes <- fm[fm$CHR == chr[i],]$genes[f] %>% str_split("___") %>% unlist()

			print(genes)
			fm[fm$CHR == chr[i],]$genes_n[f] <- length(leftout[leftout$hgnc_symbol %in% genes,]$ensembl_gene_id)
			fm[fm$CHR == chr[i],]$genes_pn[f] <- length(intersect(fmgenes, leftout[leftout$hgnc_symbol %in% genes,]$ensembl_gene_id))
			fm[fm$CHR == chr[i],]$genes_p[f] <- paste(intersect(fmgenes, leftout[leftout$hgnc_symbol %in% genes,]$ensembl_gene_id), collapse = "___")
			fm[fm$CHR == chr[i],]$genes_total[f] <- length(fmgenes)	

		}

	}

	write.table(fm, paste(codes[x], "binary_proximity_coloc_adj_fine_mapping_analysis.txt", sep = "_"))
}


results <- data.frame(phen = codes, loci = NA, ploci = NA, fm_cloci = NA, fm_loci = NA)
for (x in 1:length(codes)) {

	t <- read.table(paste(codes[x], "binary_proximity_coloc_adj_fine_mapping_analysis.txt", sep = "_"))
	results[x,2] <- nrow(t)
	t <- t[t$genes_n > 0,]
	results[x,3] <- nrow(t)
	results[x,4] <- nrow(t[t$genes_pn > 0,])

	t <- read.table(paste(codes[x], "binary_proximity_fine_mapping_analysis.txt", sep = "_"))
	t <- t[t$genes_n > 0,]
	results[x,5] <- nrow(t[t$genes_pn > 0,])

}

write.table(results, "fine_mapping_overall_summary.txt")

results$left <- results$ploci - results$fm_cloci
results2 <- results[,c(1,4,6)]
results2 <- melt(results2, id.vars = "phen")
results2$variable <- factor(results2$variable, levels = c("left", "fm_cloci"))


p <- ggplot(results2, aes(x= phen, y=value, fill=variable)) + 
    geom_bar(position="stack", stat="identity", color = "black", width = 0.75)+
    scale_fill_manual(values = c("#848585", "#98D4C2"))+
    theme_bw()

pdf("fine_mapping_summary_stack.pdf", width = 6, height = 4)
print(p)
dev.off()



