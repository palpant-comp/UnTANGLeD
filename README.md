# UnTANGLeD

![GitHub release (by tag)](https://img.shields.io/github/downloads/palpant-comp/UnTANGLeD/v1.0.0/total)

# Data and code for manuscript 'Organisation of Gene Programs Revealed by Unsupervised Analysis of Diverse Gene-Trait Associations'.
   D. Mizikovsky, M. Naval Sanchez, C. Nefzger, G. Cuellar Partida, N. J. Palpant

KEY FILES: 

1. Full_MultiXcan_matrix.txt.gz - Unprocessed data containing all MultiXcan gene-trait associations 
2. Gene_Consensus_matrix.txt.gz - Consensus matrix generated according to the optimised pipeline 
3. Gene_Correlation_matrix.txt.gz - Correlation between all genes based on MultiXcan gene-trait associations

CODE DESCRIPTION:

1. Preprocessing_clustering.R - Pre-Processing, Normalisation and Initial Clustering using Seurat 
2. Consensus_matrix_calculation.R - Consensus Matrix generation
3. Hierarchical_clustering_GO_enrichment.R - Consensus Matrix compilation, hierarchical clustering and GO enrichment of clusters
4. prediction_GWAS_genes.R - Prediction of trait associated genes 
5. fine_mapping_POPS_proximity.R - Prediction of causal genes analysis
