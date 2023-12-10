library(optparse)
library(tidyverse)
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
    make_option(c("-m", "--matching_out_path"), type="character", default=NULL, 
                help="matching_out_path"),
    make_option(c("-v", "--volumncluster_path"), type="character", default=NULL, 
                help="volumncluster_path")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# parameters
piece_use = opt$piece_use

analysis_folder = opt$analysis_folder
out_path =str_glue("{analysis_folder}/out/{piece_use}")
dir.create(out_path, recursive = T, showWarnings = F)


######## 0. Loading data and info ########################################################
# Load tracking sheet 
gs4_deauth()
google_sheet_url = ""
sheet_multi = read_sheet(google_sheet_url, sheet = 'Sample_multisections') %>% 
    filter(!is.na(PASTEFileName)) %>%
    filter(PASTEpieceName == piece_use)

# Quit if empty
if(nrow(sheet_multi) == 0){
    message(str_glue('No PASTE2 data for {piece_use}'))
    message("Quit this run.")
    quit()
}

# Load seurat objects 
obj_list = pmap(sheet_multi, ~with(list(...), {
    message("loading: ", LibraryName,' ' , piece_id, ' ' , PASTEFileName)
    readRDS(SeuratObject)
})) %>% setNames(sheet_multi$LibraryName)

# Load inidvidual PASTE2 coordinates
paste2_df_list = pmap(sheet_multi, ~with(list(...), {
    message(LibraryName,' ' , piece_id, ' ' , PASTEFileName)
    read_csv(PASTE2XinHaoCoordinate) %>% 
        mutate(
            PASTEpieceName = PASTEpieceName,
            PASTEFileName = PASTEFileName,
            LibraryName = LibraryName
        ) %>% 
    # Rename x and y into PASTE_x, PASTE_y
    dplyr::rename(PASTE_x = x, PASTE_y = y)
})) %>% setNames(sheet_multi$LibraryName) 

# Load matching spots
matching_out_path=opt$matching_out_path
if(!file.exists(str_glue('{matching_out_path}/{piece_use}_min_dist.csv'))){
    message(str_glue('{matching_out_path}/{piece_use}_matching_spot.tsv not exist'))
    message("Quit this run.")
    quit()
}
matching_df = read_csv(str_glue('{matching_out_path}/{piece_use}_min_dist.csv'))
# Add keep those names with "_2" in there
# turn into list, split by LibraryName_1 (library name of the upper slice)
matching_df_list = matching_df %>% split(f = .$LibraryName_1) %>%
    map(~as.data.frame(.) %>% 
        column_to_rownames('barcode_1') %>% 
        select(ends_with('_2')) %>% 
        setNames(str_replace(colnames(.), '_2','_NextSlice')) %>% 
        dplyr::rename(
            PASTE_x_NextSlice = x_NextSlice,
            PASTE_y_NextSlice = y_NextSlice
        )
        )

## Load Morph output
morph_df = pmap(sheet_multi, ~with(list(...), {
        read_tsv(MorphOutput) %>% 
        dplyr::rename(barcode = '...1') %>%
        # Clean columns
        select(barcode, `Filtered tumor regions`) %>%
        dplyr::rename(Filtered_tumor_regions= `Filtered tumor regions`) %>%
        mutate(
            PASTEFileName = PASTEFileName,
            PASTEpieceName = PASTEpieceName,
            LibraryName = LibraryName
            )
    })) %>% bind_rows() 


## LOAD 3D cluster result!
volumn_cluster_path = opt$volumncluster_path
if(!file.exists(str_glue('{volumn_cluster_path}/{piece_use}/5_3D_cluster_label.tsv'))){
    message(str_glue('{volumn_cluster_path}/{piece_use}/5_3D_cluster_label.tsv not exist'))
    message("Quit this run.")
    quit()
}
volumn_cluster_df = read_tsv(str_glue('{volumn_cluster_path}/{piece_use}/5_3D_cluster_label.tsv')) %>% distinct()

morph_3d_cluster_df = left_join(morph_df, volumn_cluster_df, by = c('LibraryName', 'Filtered_tumor_regions', 'PASTEFileName'))

################### Analysis ################################################################
# 1. Add meta data individual sample and save the meta.dat
# 1a. Add PASTEpieceName, PASTEFileName, PASTE_x, PASTE_y using paste2_df_list
obj_list_w_meta = obj_list
for(library_use in names(obj_list)){
    message(library_use)
    # 1a Add PASTEpieceName, PASTEFileName, PASTE_x, PASTE_y using paste2_df_list
    paste_df_use = paste2_df_list[[library_use]] %>% as.data.frame %>% column_to_rownames('barcode')
    obj_list_w_meta[[library_use]] = AddMetaData(obj_list_w_meta[[library_use]], metadata = paste_df_use)
    
    # 1b. Add 3D column cluster uinsg volumn_cluster_df
    volumn_cluster_df_use = morph_3d_cluster_df %>% filter(LibraryName == library_use) %>% 
        select(barcode, Filtered_tumor_regions, volumn_cluster) %>%
        as.data.frame %>% column_to_rownames('barcode')
    obj_list_w_meta[[library_use]] = AddMetaData(obj_list_w_meta[[library_use]], metadata = volumn_cluster_df_use)

    # 1c. Add matching barcode from adjacent section using matching_df_list
    #!! NOTE. Will not have info for the last slice!!
    #skip if library_use not in the matching_df_list
    if(!library_use %in% names(matching_df_list)){
        message(str_glue('{library_use} not in matching_df_list'))
        next
    }
    matching_df_use = matching_df_list[[library_use]]
    obj_list_w_meta[[library_use]] = AddMetaData(obj_list_w_meta[[library_use]], metadata = matching_df_use)
    obj_list_w_meta[[library_use]] = obj_list_w_meta[[library_use]]
}

# save meta.data only
iwalk(obj_list_w_meta, function(obj, library_use){
    message(library_use)
    metadata = obj@meta.data %>% rownames_to_column('barcode_orig')
    write_tsv(metadata, str_glue('{out_path}/1_PASTE_3Dcluster_metadata_{library_use}.tsv'))

})

# 2A. Creat merged object - ALL
# First check if file already exists. skip if exists
if(file.exists(str_glue('{out_path}/2_merged_obj.qs'))){
    message(str_glue('{out_path}/2_merged_obj.qs already exists'))
    message("Skip the 2A whole object merge.")
}else{
    message("merging objects")
    obj_merged = merge(
        obj_list_w_meta[[1]], obj_list_w_meta[-1], 
        add.cell.ids = names(obj_list_w_meta),
        project = piece_use
        )
    # rerun analysis
    message("rerun analysis")
    DefaultAssay(obj_merged) <- "SCT"
    VariableFeatures(obj_merged) <- map(obj_list_w_meta, VariableFeatures) %>% reduce(union)
    obj_merged <- RunPCA(obj_merged, verbose = FALSE)
    obj_merged <- FindNeighbors(obj_merged, dims = 1:30)
    obj_merged <- FindClusters(obj_merged, verbose = FALSE)
    obj_merged <- RunUMAP(obj_merged, dims = 1:30)

    # Post processing after merge: 
    # reorder libraryName. Sorted library Name in the 'names' of the vector lib_names_sorted
    lib_names_sorted = unique(obj_merged@meta.data$LibraryName) %>% 
        setNames(str_extract(.,'U[0-9]+'),.) %>% 
        gtools::mixedsort()
    # Fix NA in volumn_cluster
    obj_merged@meta.data = obj_merged@meta.data %>% 
        mutate(volumn_cluster = ifelse(is.na(volumn_cluster), 'NA', volumn_cluster))

    # Save object 
    qs::qsave(obj_merged, str_glue('{out_path}/2_merged_obj.qs'), nthreads = 20)
}

# 2B. Creat merged object - Tumor Only
if(file.exists(str_glue('{out_path}/2_merged_obj_tumor.qs'))){
    message(str_glue('{out_path}/2_merged_obj_tumor.qs already exists'))
    message("Skip the 2B tumor object merge.")
}else{
    message("merging TUMOR objects")
    obj_list_tumor = obj_list_w_meta %>% imap(function(obj, sample_id){
        tumor_cells = obj@meta.data %>% filter(Filtered_tumor_regions != '0') %>% rownames()
        print(sample_id);print(length(tumor_cells))
        # Note if morph failed, tumor_cells count will be 0
        # Return NULL to ignore the sample
        if(length(tumor_cells) == 0){ return(NULL) }
        obj = subset(obj, cells = tumor_cells) %>% FindVariableFeatures()
        return(obj)
    }) %>% compact()
    print(str_c("Samples with tumor: ", toString(names(obj_list_tumor))))
    print(str_c("Samples without tumor: ", toString(setdiff(names(obj_list_w_meta), names(obj_list_tumor)))))

    message("merging TUMOR objects")
    obj_tumor = merge(
        obj_list_tumor[[1]], obj_list_tumor[-1], 
        add.cell.ids = names(obj_list_tumor),
        project = piece_use
        )
    # rerun analysis
    message("rerun analysis")
    DefaultAssay(obj_tumor) <- "SCT"
    VariableFeatures(obj_tumor) <- map(obj_list_tumor, VariableFeatures) %>% reduce(union)
    obj_tumor <- RunPCA(obj_tumor, verbose = FALSE)
    obj_tumor <- FindNeighbors(obj_tumor, dims = 1:30)
    obj_tumor <- FindClusters(obj_tumor, verbose = FALSE)
    obj_tumor <- RunUMAP(obj_tumor, dims = 1:30)

    # Post processing after merge: 
    # reorder libraryName. Sorted library Name in the 'names' of the vector lib_names_sorted
    lib_names_sorted = unique(obj_tumor@meta.data$LibraryName) %>% 
        setNames(str_extract(.,'U[0-9]+'),.) %>% 
        gtools::mixedsort()
    # Fix NA in volumn_cluster
    obj_tumor@meta.data = obj_tumor@meta.data %>% 
        mutate(volumn_cluster = ifelse(is.na(volumn_cluster), 'NA', volumn_cluster))

    message("saving Tumor object")
    qs::qsave(obj_tumor, str_glue('{out_path}/2_merged_obj_tumor.qs'), nthreads = 20)
    }