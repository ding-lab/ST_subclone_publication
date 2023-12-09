# conda activate ROGUE
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
library(Seurat)
library(googlesheets4)

# Parameters
PATH_TO_ANALYSIS = ""
sample_info_google_sheet_url = ""
# load sheet
gs4_deauth()
sheet = read_sheet(sample_info_google_sheet_url)
sheet_use = sheet %>% filter(!is.na(AnnotationFinal))

# Load all data
st_list = pmap(sheet_use, function(LibraryName, SeuratObject, AnnotationFinal,...){
    message("Loading ", LibraryName)
    readRDS(SeuratObject)
}) %>% setNames(sheet_use$LibraryName) 

# Load and add microregion annotation
st_list = imap(st_list, function(obj, library_name){
    # load annotation_df
    annotation_df = read_csv(sheet_use %>% filter(LibraryName == library_name) %>% pull(AnnotationFinal) )
    message("Adding annotation to ", library_name)
    obj = annotation_df %>% 
        column_to_rownames('Barcode') %>% 
        .[colnames(obj),,drop=F] %>%
        setNames('Microregion') %>% 
        {AddMetaData(obj, metadata = .)}
    return(obj)
})

# subset to tumor only
tumor_list = st_list %>% imap(function(obj, library_name){
    message("Subsetting to tumor cells for ", library_name)
    obj = subset(obj, cells = obj@meta.data %>% .[!is.na(.$Microregion),] %>% rownames )
    return(obj)
})

### ROGUE workflow
# === function =================================================================
PATH_TO_ANALYSIS = ""
source(str_glue('{PATH_TO_ANALYSIS}/script/src/getROGUE.r'))

# RUN ROGUE for each sample ---------------------------------------------------------------
out_root = str_glue('{PATH_TO_ANALYSIS}/out/')

ptm = Sys.time()
rogue_list = imap(tumor_list_use, function(obj, library_name){
    # first check if output exists, if so skip
    sample_out_path = str_glue('{out_root}/{library_name}/1_ROGUE_result.rds')
    if(file.exists(sample_out_path)){
        message("ROGUE result exists for ", library_name, ". Skip the calculation.")
        return(readRDS(sample_out_path))
    }
    message("Calculating ROGUE for ", library_name)
    rogue = getROGUE(obj, assay = 'Spatial', celltype_col = 'Microregion', span = 0.6)
    # save to file
    message("Saving ROGUE for ", library_name)
    # Create folder
    sample_out_dir = str_glue('{out_root}/{library_name}/')
    dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rogue, paste0(sample_out_dir, '/1_ROGUE_result.rds'))
    return(rogue)
})

# report time in minutes ---------------------------------------------------------------
message('Total time in minutes: ', round(as.numeric(difftime(Sys.time(), ptm, units = 'mins')), 2))

# plot ROGUE score for each sample ---------------------------------------------------------------
case_col = c(
    c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"), # Accent
    c("#7FC97F", "#BEAED4", "#FDC086",  "#386CB0", "#F0027F", "#666666"), # Dark2
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF") # Set1
)

# Integrate all sample ---------------------------------------------------------------
sample_level_rouge = imap(rogue_list, function(result, name) 
	data.frame(Sample = name, ROUGE_Sample = result[['rogue.value']])) %>% 
    	bind_rows() %>% 
    	left_join(
			sheet_use[,c('SampleID_Clean','case_id','LibraryName','cancer_type','piece_section_id')], 
			by = c('Sample' = 'LibraryName')
			) 

# Add --------------------------------------------------------------------------------
# reorder case by ROUGE sample value
case_order = sample_level_rouge %>% 
	group_by(case_id) %>% 
	summarize(ROUGE_Sample = mean(ROUGE_Sample)) %>% 
	arrange(ROUGE_Sample) %>% pull(case_id) 
sample_level_rouge = sample_level_rouge %>% 
	mutate(case_id = factor(case_id, levels = case_order))

case_use = sample_level_rouge[['case_id']] %>% unique
case_col_use = case_col[1:length(case_use)] %>% setNames(case_use)

# Make dotplot -----------------------------------------------------------------------
p_rogue_cancertype = sample_level_rouge %>% 
    ggplot(aes(x = case_id, y = ROUGE_Sample, color = case_id, fill = case_id)) + 
    facet_grid(.~cancer_type, scales = 'free', space = 'free') +
    geom_violin(width = 0.5, alpha = 0.4, color = 'transparent') +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.9, aes(group = case_id),
        position = position_dodge(width = 0.5), 
        color = 'transparent'
        ) + 
    # label points with Sample name
    ggrepel::geom_text_repel(aes(label = piece_section_id), 
		size = 2, 
		position = position_dodge(width = 0.5)) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(size = 12, colour = "black", 
		angle = 90, hjust = 1, vjust = 0.5), 
          axis.title = element_text(size = 13, colour = "black")) + 
    labs(x = "Case ID", y = "ROGUE",
		title = 'Sample level heterogenity by ROGUE') + 
    # supply color
    scale_color_manual(values = case_col_use) + 
    scale_fill_manual(values = case_col_use) + 
    # remove legend
    theme(legend.position = "none")

pdf(paste0(out_root, '/1_ROGUE_cancerType.pdf'), w = 7, h = 5)
    p_rogue_cancertype %>% print()
dev.off()

