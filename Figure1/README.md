## Figure 1. Spatial microregions across cancers
### Directory structure
* `0_Load_All_Sample` : Load Visium Data cohort and add metadata
	* `/table` :  Table contained the genetic clones determined in section `Figure 2. Genetic alteration reveals spatial clonal evolution`
* `1_Calculate_metrics` : Calculate the area size per spot, split tumor microregion in to 3 different size groups, estimate tumor micoregion density per sample
* For tumor micoregion determination using `Morph` tool, plase refer to this repo: https://github.com/ding-lab/morph