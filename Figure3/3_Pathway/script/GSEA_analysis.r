# conda activate clusterprofiler
# 2023/09/13 Simon Mo
### ------------------ Load data ------------------ ###
# Load samples using script
source('/PATH/TO/0_Load_All_Sample/script/load_samples.r')

# Parameters
out_path = ''

### ------------------ Function ------------------ ###
src_path="/PATH/TO/3_Pathway/script"
source(str_glue('{str_path}/function_findmarker_enhanced.r'))
# Load GSEA pathway analysis scripts
source(str_glue('{src_path}/function_runGSEA.r'))
# Heatmap function
source(str_glue('{src_path}/function_AnnDotPlot.r'))
# Module score
source(str_glue('{src_path}/function_TestAndFixAddModuleScore.r'))
# Add average expression of genetic subclone 
source(str_glue('{src_path}/function_CalculateCloneExpressionLevel.r'))
# Load data base
LoadMSigDBHuman()

### ------------------ Analysis ------------------ ###
### ------------------ Find DEG, GSEA ------------------ ###
# plot 
library(Seurat)
library(enrichplot)
library(patchwork)

out_path_analysis = str_glue('{out_path}/1A_subclone_vs_TME_plot')
#ST_microregion_deg_list
st_groupby = 'genetic_clone'

## -------- This GSEA workflow is generic and can be apply to many datasets -------- ##
# Get DEGs, GSEA for each sample
sample_use_all = names(st_list)
for(sample_use in sample_use_all){
    # 0. select sample
    message("Processing sample:", sample_use)
    ST_use = st_list[[sample_use]]
    Idents(ST_use) = st_groupby
    #ST_use = subset(ST_use, downsample = 5) # For Testing

    # 1. Run DEG analysis
    # TME as ref, run A. miroregion vs TME, B. TME vs NotTME
    deg_all_df = FindMarkersEachVsRefComplete(ST_use, group.by = st_groupby, ident_ref = '0')
    
    # save result
    dir.create(str_glue('{out_path_analysis}/{sample_use}'), recursive = TRUE, showWarnings = FALSE)
    write_tsv(deg_all_df, str_glue('{out_path_analysis}/{sample_use}/0_DEG_result.tsv'))

    # 2. get GSEA result for TME vs Tumor for multiple genesets
    gsea_genesets_list = FindAllMarkerTable2GSEAresult_MsigDB(deg_all_df , genesets_include = c("H","C6")) %>%
        discard(.p = ~length(.)==0)

    # 3. plot GSEA result
    iwalk(gsea_genesets_list, function(result_list, geneset_name){
        # save result
        dir.create(str_glue('{out_path_analysis}/{sample_use}'), recursive = TRUE, showWarnings = FALSE)
        saveRDS(result_list, str_glue('{out_path_analysis}/{sample_use}/1_GSEA_result_{geneset_name}.rds'))
        # plot
        p_st = SpatialDimPlot(ST_use, group.by = st_groupby, stroke = NA, image.alpha = 0, label = T)
        pdf(str_glue('{out_path_analysis}/{sample_use}/1_dotplot_{geneset_name}.pdf'), height = 12, width = 10)
            MakeGeneSetDotplot(result_list) %>% print()
            print(p_st)
        dev.off()
    })
}

# ------------------ Set up parameters ------------------ #
markerset_name = "H"

# ------------------ PLOT ------------------ #
# 1. Plot DEG heatmap
# B. Extract top GSEA pathways for nonTME and plot top features
sample_use_all = names(st_list)

n_genes_plot = Inf # set Inf to plot all genes
for(sample_use_name in sample_use_all){
    # Parameters    
    gsea_file_path = str_glue('{out_path_analysis}/{sample_use}/1_GSEA_result_{geneset_name}.rds')
    if(!file.exists(gsea_file_path)) next
    message("Processing sample:", sample_use_name)
    gsea_use_list = readRDS(gsea_file_path)
    st_obj_use = st_list[[sample_use_name]]

    # Get GSEA core genes for a list of GSEA results
    # Select 15 genes from each select GSEA result to plot
    gsea_genes_df = imap(gsea_use_list, function(gsea_use, ident){
        GetGSEAgenes(gsea_use) %>% mutate(ident = ident)
    }) %>% bind_rows() %>% 
        distinct(ID, core_enrichment) %>% 
        group_by(ID) %>% 
        slice_head(n=n_genes_plot)
    gsea_genes_plt_list = split(gsea_genes_df$core_enrichment, gsea_genes_df$ID)

    # PreCheck if will cause issue when running AddModuleScore
    st_obj_use = TestAndFixSeuratForAddModuleScore(st_obj_use)
    
    # Annotation heatmap
    dir.create(str_glue('{out_path_analysis}/{sample_use_name}/2_ExpHeatmap/{markerset_name}/'), recursive = TRUE, showWarnings = FALSE)
    p_list = imap(gsea_genes_plt_list, possibly(function(features_plt, geneset_id){
        AnnoDotPlot(st_obj_use, group.by = 'Filtered_tumor_regions', features = features_plt, 
            annotation_idents = c('genetic_clone'), label_ident= T,
            cluster_row = T,
            cluster_col = T,
            mode = 'Heatmap',
            title = geneset_id,
            subtitle = sample_use_name,
            ModuleScoreHeight = 3,
            highlight_tiles = F,
            highlight_cutoff_quantile = 0.6,
            highlight_color = "#333333",
            highlight_thickness = 0.5,
            ) 

    }, otherwise = NULL)) %>% discard(.p = ~is.null(.x)) 

    # Plot
    iwalk(p_list, function(p, geneset_id){
        message("Plotting:", geneset_id)
        pdf(str_glue('{out_path_analysis}/{sample_use_name}/2_ExpHeatmap/{markerset_name}/2_annoDotplot_{geneset_id}.pdf'), height = 8, width = 8)
            print(p)
        dev.off()
    })
}

## -----    Plot tumor region expression plot for each pathway    ----- ##
## Next Plot tumor region expression plot for each pathway
# Version 2 - 20231010

# Extract the genes list and save
# For loop and split gene by panel
samples_use = names(st_list)

for(sample_use_name in samples_use){
    # Parameters    
    gsea_file_path = str_glue('{out_path_analysis}/{sample_use}/1_GSEA_result_{geneset_name}.rds')
    if(!file.exists(gsea_file_path)) next
    message("Processing sample:", sample_use_name)
    gsea_use_list = readRDS(gsea_file_path)
    st_obj_use = st_list[[sample_use_name]]
    st_tumor_use = tumor_list[[sample_use_name]]
    
    panels_per_file = 9
    gsea_genes_df = imap(gsea_use_list, function(gsea_use, ident){
        GetGSEAgenes(gsea_use) %>% mutate(ident = ident)
    }) %>% bind_rows() %>% 
        #filter(ID %in% gsea_id_plt) %>%  # Use all geneset
        distinct(ID, core_enrichment) %>% 
        group_by(ID) %>% 
        # Split by number of panels
        mutate(ID_split = ceiling(seq_along(core_enrichment)/panels_per_file)) %>% 
        mutate(ID_split_full = str_c(ID, '_', ID_split)) 
    
    # Add average expression of genetic subclone 
    gsea_genes_df = gsea_genes_df %>% 
    left_join(
        y = CalculateCloneExpressionLevel(obj_use, features = unique(.$core_enrichment)),
        by = c('core_enrichment' = 'Gene')
    ) %>% # rearrange by clone group
        arrange(ID, max_tumor)
    write_tsv(gsea_genes_df, str_glue('{out_path_analysis}/{sample_use_name}/4_GSEA_result_long_{markerset_name}.tsv'))
    
    # filtered version 
    gsea_genes_filtered_df = gsea_genes_df %>% filter(min_max_tumor_ratio > 1.5, max_tme_ratio > 1.5) %>%
        # Rearrange
        arrange(ID, max_tumor) %>%
        # Redo splitting 
        group_by(ID) %>% 
        # Split by number of panels
        mutate(ID_split = ceiling(seq_along(core_enrichment)/panels_per_file)) %>% 
        mutate(ID_split_full = str_c(ID, '_', ID_split)) %>%
        mutate(ID_split_full = str_c(ID_split_full, '_', max_tumor)) 

    write_tsv(gsea_genes_filtered_df, str_glue('{out_path_analysis}/{sample_use_name}/4_GSEA_result_long_{markerset_name}_filtered.tsv'))


    # Plot
    gsea_genes_plt_list = split(gsea_genes_filtered_df$core_enrichment, gsea_genes_filtered_df$ID_split_full)
    # SpatialPlot
    dir.create(str_glue('{out_path_analysis}/{sample_use_name}/5_SpatialPltPathwayGenesFiltered/{markerset_name}/'), recursive = TRUE, showWarnings = FALSE)
    iwalk(gsea_genes_plt_list, function(features_plt, geneset_name){
        p = SpatialPlot(st_tumor_use, features = features_plt, stroke = NA, image.alpha = 0.4)
        pwhole = SpatialPlot(st_obj_use, features = features_plt, stroke = NA, image.alpha = 0.4)
        message("Plotting:", geneset_name)
        pdf(str_glue('{out_path_analysis}/{sample_use_name}/5_SpatialPltPathwayGenesFiltered/{markerset_name}/3_SpatialFeature_{geneset_name}.pdf'), height = 8, width = 8)
            print(p)
            print(pwhole)
        dev.off()
    })
    }
