# conda activate seurat4.3
library(tidyverse)
library(Seurat)
library(googlesheets4)

# Load subclone assignment sheet 
gs4_deauth()
# 2023/09/07: update to include more sample
genetic_clone_sheet = read_tsv('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_16-InferCNV_OCT_workflow//OCT_genetic_clones_2023-09-07.tsv')

# Load sample sheet  
sample_sheet_url = ""
sample_sheet = read_sheet(, sheet = 'Sample_level')
sheet_use = sample_sheet %>% filter(LibraryName %in% genetic_clone_sheet$sample_id)

# Load all data
st_list = pmap(sheet_use, function(LibraryName, SeuratObject, MorphOutput,...){
    message("Loading ", LibraryName)
    readRDS(SeuratObject)
}) %>% setNames(sheet_use$LibraryName) 

# Add annotation
st_list = imap(st_list, function(obj, library_name){
    # load annotation_df
    annotation_df = read_tsv(sheet_use %>% filter(LibraryName == library_name) %>% pull(MorphOutput) ) %>% 
        as.data.frame() %>% 
    # move first column to the barcode
        column_to_rownames('...1') %>% 
        .[,c('Filtered tumor regions'), drop =F]
    message("Adding annotation to ", library_name)
    obj = annotation_df %>% 
        .[colnames(obj),,drop=F] %>%
        setNames('Filtered_tumor_regions') %>% 
        {AddMetaData(obj, metadata = .)}
    return(obj)
})

# Add genetic clone 
st_list = imap(st_list, function(obj, sample_ID){
    genetic_sheet_use = filter(genetic_clone_sheet, sample_id == sample_ID)
    obj@meta.data$genetic_clone = plyr::mapvalues(
        x = obj@meta.data$Filtered_tumor_regions,
        from = genetic_sheet_use$Filtered_tumor_regions,
        to = genetic_sheet_use$genetic_clone
    )
    return(obj)
})

# subset to tumor only
tumor_list = st_list %>% imap(function(obj, library_name){
    message("Subsetting to tumor cells for ", library_name)
    obj = subset(obj, cells = obj@meta.data %>% filter(!Filtered_tumor_regions == 0) %>% rownames )
    return(obj)
})

# Report object that can be use
message("Avaialbe objects: ")
message("tumor only list: tumor_list")
message("all cells list: st_list")

