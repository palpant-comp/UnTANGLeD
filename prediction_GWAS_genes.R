### GWAS PREDICTION ###


library(cluster)
library(stringr)
library(ggplot2)

  full_data <- "/home/uqdmizik/90days/gwas_test/pipeline/"
  multixcan <- "/home/uqdmizik/90days/gwas_test/multixcan/REAL_PHEN_OUTPUTS/multixcan/"

  tog <- readRDS(paste0(full_data, "full_matrix.rds"))

  phenotypes <- list.files(multixcan)
  phenotypes <- phenotypes[grep("old|mid|asian", phenotypes)]
  phenotypes[7] <- "SCZ_new_multixcan.txt"
  pcodes <- c("gout", "hdlasian", "hdl", "ldlasian", "ldl", "ra", "scz_asian", "sczold", "sle_asian", "tgasian", "tg", "uc", "base")
  newpcodes <- c(1, 4, 4, 8, 8, 10, 12, 12, 15, 18, 18, 20, 6)
  inflection = 242
  mcomp <- matrix(nrow = length(pcodes), ncol = 10)

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


  for (x in 1:length(pcodes)) {

	t <- tog
	t3 <- as.matrix(t)
	t3[is.na(t3)] <- 1


	t3 <- t3[,colnames(t3) %in% colnames(cm)]

	t3 <- -log10(t3)
	t3[is.infinite(t3)] <- 315

	if(x == 13) { } else {

	mxp <- read.table(paste0(multixcan, phenotypes[x]), head = T)
	mxp$gene <- str_split_fixed(mxp$gene, "\\.", 2)[,1]


	print("data loaded")

	mxp <- mxp[mxp$gene %in% colnames(t),]
	mxp2 <- mxp$pvalue
	names(mxp2) <- mxp$gene

	adt <- rep(1, length(setdiff(colnames(t), names(mxp2))))
	names(adt) <- setdiff(colnames(t), names(mxp2))
	mxp2 <- c(mxp2, adt)
	mxp2 <- mxp2[order(match(names(mxp2),colnames(t)))]
	t <- rbind(t, mxp2)
	rownames(t)[1394] <- "NEW_PHENOTYPE"
	print("data integrated")

}

	t3 <- as.matrix(t)
	t3[is.na(t3)] <- 1


	t3 <- t3[,colnames(t3) %in% colnames(cm)]

	t3 <- -log10(t3)
	t3[is.infinite(t3)] <- 315

	cluster_stats <- read.table("/scratch/90days/uqdmizik/gwas_test/globalcc/old/clus242_metrics_summary.txt")
	cluster_stats <- cluster_stats[,1:3]


		ggplot(cluster_stats, aes(x = silhouette, y = correlation))+geom_point()



	if(x == 13) {

			npm <- mean(t3[1034,])

	    np <- t3[1034,]




	} else {

		npm <- mean(t3[1394,])

	  np <- t3[1394,]
	}

	for (i in 1:inflection) {
  		cluster <- names(sub_grp[sub_grp == i])
  		cluster_stats[i,4] <- mean(np[names(np) %in% cluster])/npm
  		cluster_stats[i,5] <- t.test(np[names(np) %in% cluster], np[!(names(np) %in% cluster)])$p.value
	}

	files <- list.files(multixcan)
	newg <- read.table(paste0(multixcan, files[newpcodes[x]]), header = T)
	newg$gene <- str_split_fixed(newg$gene, "\\.", 2)[,1]
	nsigg <- na.omit(newg[newg$pvalue < 1e-4,]$gene)
	nsigg <- nsigg[nsigg %in% colnames(t3)]



	cluster_stats$nps <- NA
	cluster_stats$npps <- NA

	cluster_stats$nps <- ifelse(cluster_stats[,4] > 1.5, "S", "NS")

	for (i in 1:inflection) {
  
  		if (cluster_stats[i,5] < 0.05 && cluster_stats[i,4] > 1) {
    		cluster_stats$npps[i] <- "S"
  		} else {
    		cluster_stats$npps[i] <- "NS"
  		}
  
	}


if(x == 13) {

	osig <- names(np[np > 4])

	} else {
	osig <- na.omit(mxp[mxp$pvalue < 1e-4,]$gene)

	}


colnames(cluster_stats)[1:5] <- c("cluster", "sil", "cor", "mean_p", "p_sig")
cluster_stats$sigg <- ifelse(cluster_stats$cluster %in% sub_grp[names(sub_grp) %in% osig], "sig", "nsig")
cluster_stats$comb1 <- paste0(cluster_stats$nps, cluster_stats$npps)
p <- ggplot(cluster_stats, aes(x = sil, y = cor, color = comb1, size = mean_p))+geom_point()+theme_minimal()


asg <- names(sub_grp[sub_grp %in% cluster_stats[cluster_stats$comb1 %in% c("NSS", "SNS", "SS"),]$cluster])
op <- intersect(nsigg, asg)

	## FOUR - NPS OR/AND NPPS (EITHER) ##

	expected <- length(intersect(osig, op))
	sub_grp[names(sub_grp) %in% op]

	ctg <- matrix(nrow = 2, ncol =2)
	colnames(ctg) <- c("cluster", "ncluster")
	rownames(ctg) <- c("gwas", "ngwas")
	ctg[1,1] <- length(op) - expected
	ctg[2,1] <- length(asg) - length(op)
	ctg[1,2] <- length(intersect(nsigg, colnames(t3)[!(colnames(t3) %in% asg)]))
	ctg[2,2] <- length(colnames(t3)) - ctg[1,1] - ctg[2,1] - ctg[1,2]

	mcomp[x,1] <- ctg[1,1]
	mcomp[x,2] <- ctg[1,2]
	mcomp[x,3] <- ctg[2,1]
	mcomp[x,4] <- ctg[2,2]
	
	if (chisq.test(ctg)$expected[1,1] > 5) {
		mcomp[x,5] <- chisq.test(ctg)$p.value
		mcomp[x,10] <- "chisq"
	} else {

		mcomp[x,5] <- fisher.test(ctg)$p.value
		mcomp[x,10] <- "fisher"
	}
}

	mcomp[x,6] <- length(osig)
	mcomp[x,7] <- length(nsigg)
	mcomp[x,8] <- expected
	mcomp[x,9] <- length(intersect(osig, nsigg))


	trial <- as.data.frame(sub_grp)

	trial$gwas <- ifelse(rownames(trial) %in% nsigg, "GWAS", "NOT")
	trial$cluster <- ifelse(trial$sub_grp %in% sub_grp[names(sub_grp) %in% asg], "CLUSTER", "NOT")
	trial <- trial[!(rownames(trial) %in% osig),]
	p <- ggplot(trial, aes(x = gwas, fill = cluster))+geom_bar(position='fill')+theme_bw()

	pdf(paste(pcodes[x], "proportion_plot_gwclus.pdf", sep = "_"))
	print(p)
	dev.off()

	p <- ggplot(trial, aes(fill = gwas, x = cluster))+geom_bar(position='fill')+theme_bw()

	pdf(paste(pcodes[x], "proportion_plot_clusgw.pdf", sep = "_"))
	print(p)
	dev.off()

	trial2[x,1] <- pcodes[x]
	trial2[x,2] <- ctg[1,1]
	trial2[x,3] <- ctg[1,2]
	trial2[x,4] <- ctg[2,1]
	trial2[x,5] <- ctg[2,2]

cluster_stats$new <- NA
cluster_stats$newn <- NA

cluster_stats[cluster_stats$cluster %in% names(table(sub_grp[names(sub_grp) %in% nsigg])),]$new <- "New"

ev <- rep(0, inflection-length(table(sub_grp[names(sub_grp) %in% nsigg])))
names(ev) <- cluster_stats$cluster[!(cluster_stats$cluster %in% names(table(sub_grp[names(sub_grp) %in% nsigg])))]
ev <- c(ev, table(sub_grp[names(sub_grp) %in% nsigg]))

for (i in 1:inflection) {
  cluster_stats[i,]$newn <- ev[names(ev) == i]
}

cluster_stats$comb2 <- ifelse(cluster_stats$comb1 %in% c("NSS", "SNS", "SS"), "ass", "not_ass")
cluster_stats$color_col <- paste0(cluster_stats$comb2, cluster_stats$new)
cluster_stats[cluster_stats$newn == 1,]$newn <- 0

#p <- ggplot(cluster_stats, aes(x = sil, y = cor, size = newn, color = color_col))+geom_point()+theme_minimal()

#pdf(paste(pcodes[x], "prediction_plot.pdf", sep = "_"), width = 6.5, height = 6)
#print(p)
#dev.off()


  }


  rownames(mcomp) <- pcodes

  colnames(mcomp) <- c("ctg11", "ctg12", "ctg21", "ctg22", "enrichment", "old_genes", "new_genes", "expected", "overlap", "stat_test")

  write.table(mcomp, "enrichment_result_summary_RC_no_integration.txt")

  colnames(trial2) <- c("phen", "gw_clus", "gw_nclus", "ngw_clus", "ngw_nclus")
  t2 <- melt(trial2, id=c("phen", "gw_nclus", "gw_clus"))
  t2 <- t2[,c(1, 4:5)]
  t2$phen <- paste(t2$phen, "ngw", sep = "_")
  t3 <- melt(trial2, id=c("phen", "ngw_clus", "ngw_nclus"))
  t3 <- t3[,c(1, 4:5)]
  t3$phen <- paste(t3$phen, "gw", sep = "_")
  t2 <- rbind(t3,t2)
  t2 <- t2[order(t2$phen),]


  p <- ggplot(t2, aes(fill=variable, y=value, x=phen)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values = c("#98D4C2", "#848585", "#98D4C2", "#848585"))+
    theme_bw()

  pdf(paste(pcodes[x], "proportion_plot_clusgw.pdf", sep = "_"), width = 15, height = 5)
	print(p)
	dev.off()



  "enrichment_result_summary_mat_integration.txt"