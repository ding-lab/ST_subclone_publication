## Figure 6. 3D tumor volume reconstruction revealed diverse tumor growth patterns
### Directory structure
* `0_data` :  Table contained the genetic clones determined in section `Figure 2. Genetic alteration reveals spatial clonal evolution`
* `0_palette` : Color palette used in Sankey Diagram
* `1_run_PASTE2` : workflow and script to run PASTE2 for multisection alignment
* `2_PASTE2_analysis` : Analysis scripts after PASTE2 alignment to get spot level connections and metrics calculation for reconstructed 3D volumes
	- `1_plot_paste_result` : Visualize paste alginment result
	- `2_get_matching_spot` : Identify closest spot between adjacent slice for 3D volume reconstrunction
	- `3_multi_section_regions_N_min_3` : Make Sankey plot to visualize connections microregions among multisections. Determine 3D clusters
	- `4_mergeobj_and_3dcluster_N_min_3` : Visiualize 3D cluster using Seurat Spatial Plots
	- `5_calculate_metrics` : Calculate and visualize 3D structure related metrics: Degree, Loop, N component 
