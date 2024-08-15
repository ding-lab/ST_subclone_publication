library(optparse)
library(tidyverse)
library(patchwork)
library(googlesheets4)
library(tictoc)
library(future)
library(furrr)
library(Seurat)
library(qs)

plan(multisession, workers = 20)
options(future.globals.maxSize = Inf)

## Write a optparse list to input necessory info
option_list = list(
    make_option(c("-a", "--analysis_folder"), type="character", default=NULL, 
                help="analysis_folder"),
    make_option(c("-p", "--piece_use"), type="character", default=NULL, 
                help="piece_use"),
    make_option(c("-m", "--merged_obj_path"), type="character", default=NULL, 
                help="merged_obj_path")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# parameters
#piece_use = 'ht268'
piece_use = opt$piece_use

#analysis_folder='/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/30-PASTE2/1_PASTE2_result/5_mergeobj_and_3dcluster'
analysis_folder = opt$analysis_folder
out_path =str_glue("{analysis_folder}/out/{piece_use}")
dir.create(out_path, recursive = T, showWarnings = F)

# volumn_folder
# volumn_folder="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/30-PASTE2/1_PASTE2_result/5_mergeobj_and_3dcluster/out"
merged_folder=opt$merged_obj_path

######## 0. Loading data and info ########################################################
# Load tracking sheet 
gs4_deauth()
sheet_multi = read_sheet('https://docs.google.com/spreadsheets/d/1MhZ98AVoUQoUxIiE4pX4pbL-eEy-RQeF9kV_xLwRz9E/edit#gid=1451401776', sheet = 'Sample_multisections') %>% 
    filter(!is.na(PASTEFileName)) %>% 
    filter(PASTEpieceName == piece_use)

# Load merged object
merged_path = str_glue('{out_path}/2_merged_obj.qs')
merged_tumor_path = str_glue('{out_path}/2_merged_obj_tumor.qs')
obj_merged = qs::qread(merged_path, nthreads = 10)
obj_tumor = qs::qread(merged_tumor_path, nthreads = 10)

######## A. Make Plots ########################################################
# Starting from Number 3
n_slice = length(Images(obj_merged))

#### Preprocess
# reorder libraryName. Sorted library Name in the 'names' of the vector lib_names_sorted
lib_names_sorted = unique(obj_merged@meta.data$LibraryName) %>% 
    setNames(str_extract(.,'U[0-9]+'),.) %>% 
    gtools::mixedsort()
# Fix NA in volumn_cluster
obj_merged@meta.data = obj_merged@meta.data %>% 
    mutate(volumn_cluster = ifelse(is.na(volumn_cluster), 'NA', volumn_cluster))

#### PLOT
obj_merged@meta.data$LibraryName = factor(obj_merged@meta.data$LibraryName, 
    levels = names(lib_names_sorted))

# 3a. SpatialDimPlot
idents_plot = c('seurat_clusters', 'Filtered_tumor_regions', 'volumn_cluster')
for(ident_use in idents_plot){
    message(ident_use)
    p_st = SpatialPlot(obj_merged, group.by = ident_use, 
    stroke = NA, crop = F, image.alpha = 0.8,
    label = T, label.box = F, label.color = 'black', combine =F,
    ) %>% wrap_plots(guides = 'collect', nrow = 1) & theme(legend.position = 'bottom') & NoLegend()
    pdf(str_glue('{out_path}/3a_{ident_use}_SpatialDimPlot.pdf'), width = 20, height = 5)
    print(p_st)
    dev.off()

    p_umap_merged = DimPlot(obj_merged, group.by = ident_use, raster = T) 
    p_umap_split = DimPlot(obj_merged, group.by = ident_use, split.by = 'LibraryName', raster = T, combine =F, label = T, repel = T)
    p_umap_all = c(list(p_umap_merged), p_umap_split) %>% 
        wrap_plots(guides = 'collect', nrow = 1, widths = c(1, n_slice)) & 
        theme(legend.position = 'bottom') & coord_fixed()
    pdf(str_glue('{out_path}/3a_{ident_use}_umap.pdf'), width = 20, height = 5)
    print(p_umap_all)
    dev.off()
}


######## B. Make Plots for tumor only ########################################################
# reorder libraryName. Sorted library Name in the 'names' of the vector lib_names_sorted
lib_names_sorted = unique(obj_tumor@meta.data$LibraryName) %>% 
    setNames(str_extract(.,'U[0-9]+'),.) %>% 
    gtools::mixedsort()
# Fix NA in volumn_cluster
obj_tumor@meta.data = obj_tumor@meta.data %>% 
    mutate(volumn_cluster = ifelse(is.na(volumn_cluster), 'NA', volumn_cluster))

# 3b. SpatialDimPlot on tumor only
idents_plot = c('seurat_clusters', 'Filtered_tumor_regions', 'volumn_cluster')
#idents_plot ='seurat_clusters'
for(ident_use in idents_plot){
    message(ident_use)
    p_st = SpatialPlot(obj_tumor, group.by = ident_use, 
    stroke = NA, crop = F, image.alpha = 0.8,
    label = T, label.box = F, label.color = 'black', combine =F,
    ) %>% wrap_plots(guides = 'collect', nrow = 1) & theme(legend.position = 'bottom')
    pdf(str_glue('{out_path}/3b_tumor_{ident_use}_SpatialDimPlot.pdf'), width = 20, height = 5)
    print(p_st)
    dev.off()

    p_umap_merged = DimPlot(obj_tumor, group.by = ident_use, raster = T) 
    p_umap_split = DimPlot(obj_tumor, group.by = ident_use, split.by = 'LibraryName', raster = T, combine =F, label = T, repel = T)
    p_umap_all = c(list(p_umap_merged), p_umap_split) %>% 
        wrap_plots(guides = 'collect', nrow = 1, widths = c(1, n_slice)) & 
        theme(legend.position = 'bottom') & coord_fixed()
    pdf(str_glue('{out_path}/3b_tumor_{ident_use}_umap.pdf'), width = 20, height = 5)
    print(p_umap_all)
    dev.off()
}
