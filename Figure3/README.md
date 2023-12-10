## Figure 3. Genetic alteration drives transcriptional similarity in tumor microregions
### Directory structure
* `0_Load_All_Sample` : Load Visium Data cohort and add metadata
	* `/table` :  Table contained the genetic clones determined in section `Figure 2. Genetic alteration reveals spatial clonal evolution`
* `1_Heterogenity_score` : Calculate the heterogenity score : 1 - ROGUE
* `2_Correlation` : Calculate and visualize tumor microregion correlations
* `3_Pathway` : GSEA pathway analysis for tumor micoregions vs TME regions