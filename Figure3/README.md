## Scripts related to "Genetic changes drive tumour disparities" & "Cellular pathway at tumour core and edge" (Figure 3)
### Directory structure
* `0_Load_All_Sample` : Load Visium Data cohort and add metadata
	* `/table` :  Table contained the genetic clones determined in section `Figure 2. Genetic alteration reveals spatial clonal evolution`
* `1_Heterogenity_score` : Calculate the heterogenity score : 1 - ROGUE
* `2_Correlation` : Calculate and visualize tumor microregion correlations
* `3_Pathway` : GSEA pathway analysis for tumor micoregions vs TME regions
* `4_Layey` : Scripts related to "Cellular pathway at tumour core and edge"
** `1_layer_statistics.R`: Summarize number of layers and size in each tumor microregion across the cohort (related to Fig. 4a-h)
** `2_layer_DEG_correlation.R`: Identify gene expression correlated with tumor depths (in layers) and their associated biological pathways in each tumor sample (related to Fig. 4c-h) 
** `3_layer_DEG_pathway_analysis.R`: Cohort-level summary of genes recurrently associated with tumor center/periphery and their associated biological pathways (related to Fig. 4f-h)



