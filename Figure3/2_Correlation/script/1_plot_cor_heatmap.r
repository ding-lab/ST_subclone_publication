library(tidyverse)
library(Seurat)
library(googlesheets4)
library(future)
library(furrr)
library(tictoc)
library(qs)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(patchwork)


options(future.globals.maxSize = 1000 * 1024 ^ 2)
plan(multicore, workers = 10)

############################################################################
# Load ST data and add manual annotation
############################################################################
# Load samples using script

source('/PATH/TO/0_Load_All_Sample/script/load_samples.r')


# set output
out_path = ""
fig_out_root = ""
############################################################################

############################################################################
# Get Pairwise correlation of each subclone
############################################################################
# Get As many discrete color from RColorBrewer
distinct_color = RColorBrewer::brewer.pal.info %>% 
    filter(category == 'qual') %>% 
        rownames_to_column('palette') %>% 
        pmap(function(palette, maxcolors, ...){
    RColorBrewer::brewer.pal(n = maxcolors, palette)
}) %>% unlist()

# Plot
tumor_spatial_lognorm_list = imap(tumor_list, function(stobj, sample_ID){
    message("Processing: ", sample_ID)
    DefaultAssay(stobj) = 'Spatial'
    NormalizeData(stobj)
})

# Recalculate variable features
tumor_spatial_lognorm_list = imap(tumor_spatial_lognorm_list, function(stobj, sample_ID){
    message("Processing: ", sample_ID)
    stobj = FindVariableFeatures(stobj, selection.method = 'vst', nfeatures = 2000)
    return(stobj)
})

# Run all combinations
run_group = expand.grid(c('SCT','Spatial'), c('counts', 'data'), c(100, 500, 2000)) %>% 
    as.data.frame %>% 
    setNames(c('assay_use', 'slot_use', 'n_features_use')) %>% 
    mutate(across(where(is.factor), as.character)) %>% 
    as_tibble

# select sample list
sample_list_use = tumor_spatial_lognorm_list
pwalk(run_group, function(assay_use, slot_use, n_features_use){
    message("Running: assay = ", assay_use, ", slot_use = ", slot_use, " n feautres = ", n_features_use)
    # set output
    fig_out_path = str_glue("{fig_out_root}/NFeature_{n_features_use}/")
    dir.create(file.path(fig_out_path, assay_use), showWarnings = F, recursive = T)

    # Make correlation heatmap plot
    #samples_use_all = "HT112C1-U1_ST_Bn1"
    samples_use_all = names(sample_list_use)
    for(sample_use in samples_use_all){
        message("Plotting: ", sample_use)
        st = sample_list_use[[sample_use]]
        
        # 1. Get Average Expression of Variable Features
        avg_exp_subclone = AverageExpression(
                st, 
                features = VariableFeatures(st)[1:n_features_use], 
                group.by = 'Filtered_tumor_regions', 
                assays = assay_use, slot = slot_use
            ) %>% .[[1]]
        message("Average expression matrix dimension: ", dim(avg_exp_subclone))
        # 2. Get Correlation of each subclone
        cor_subclone = cor(avg_exp_subclone, method = 'pearson') %>% as.data.frame

        # 3. Get genetic_clone clone annotation and color
        genetic_clone_df = FetchData(st, vars = c("genetic_clone", "Filtered_tumor_regions")) %>% 
            distinct() %>% remove_rownames() %>% 
            column_to_rownames('Filtered_tumor_regions')
        genetic_clone_vec = genetic_clone_df$genetic_clone %>% setNames(rownames(genetic_clone_df))

        # Get subclone color + Annotate
        col_genetic_clone = colorspace::qualitative_hcl(length(unique(genetic_clone_vec)), "Dynamic") %>% setNames(unique(genetic_clone_vec))
        top_anno  = HeatmapAnnotation(
            GeneticClone = genetic_clone_vec[rownames(cor_subclone)], 
                col = list(GeneticClone = col_genetic_clone),
                show_legend = c(FALSE)
            )
        left_anno = rowAnnotation(
            GeneticClone = genetic_clone_vec[rownames(cor_subclone)], 
                col = list(GeneticClone = col_genetic_clone),
                show_legend = c(TRUE)
            )
        # 4. Plot Correlation heatmap
        text_size = 15
        cor_subclone_mtx = cor_subclone %>% as.matrix
        phm = Heatmap(cor_subclone_mtx, 
            top_annotation = top_anno,
            left_annotation = left_anno,
            name = 'Pearson\nCorrelation', 
            cluster_rows = T, cluster_columns = T,
            width = unit(15, 'cm'), height = unit(15, 'cm'),
            col = rev(RColorBrewer::brewer.pal(9, "RdYlGn")),
            column_names_gp = grid::gpar(fontsize = text_size),
            row_names_gp = grid::gpar(fontsize = text_size))

        #  7 Save plots
        dir.create(file.path(fig_out_path, assay_use, slot_use), showWarnings = F, recursive = T)
        pdf(file.path(fig_out_path, assay_use, slot_use, paste0('Cor_', sample_use, '.pdf')), width = 10, height = 8)
            print(phm)
        dev.off()
    }
    

})
