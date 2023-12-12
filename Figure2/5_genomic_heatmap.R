# Plot genomic heatmap with spatial CNV 
# Jingxian Clara Liu
# 2023/10/20

library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
library(colorspace)
library(circlize)
library(optparse)
library(infercnv, lib.loc = "/diskmnt/Projects/Users/cliu/software/anaconda3/envs/inferCNV/lib/R/library")
library(doParallel)
library(foreach)
library(EnrichedHeatmap)
library(patchwork)

set.seed(1234)

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")

# Test parameters
setting = "w151_dmeanvar_n200"
case_id = "HT112C1"
MAX_K=10
linkage_method = "ward.D2"
output_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_16-InferCNV_OCT_workflow/"

hatchet_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/Cong/hatchet_wes_v6/"

manual_K = c("HT260C1"=2, "HT397B1"=3, "HT268B1"=2, "HT112C1"=1, "HT306P1"=2, "HT270P1"=4)

### Parameters ------------------------------------------------------
option_list = list(
	  make_option(c("-c", "--case_id"),
		type="character",
		default=NULL,
		help="Sample name for ST",
		metavar="character"),
	  make_option(c("-p", "--setting"),
		type="character",
		default=setting,
		help="Setting parameter for InferCNV",
		metavar="character")
	  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

case_id=opt$case_id
setting=opt$setting

### Set up sample to be processed -----------------------------------
selected_tbl = read_tracking_sheet()
sample_id_list = selected_tbl %>% filter(Study_cohort=="Discovery", case_id==!!case_id) %>% pull(Sample) %>% sort
WES_sample_id = read_column_val(sample_id_list[1], "WES_sample_id")
calico_sample_id = read_column_val(sample_id_list[1], "calico_sample_id")
sn_sample_id = read_column_val(sample_id_list[1], "Paired_snRNA_id")
hatchet_run = grep(case_id, dir(hatchet_dir), value=T)
stopifnot(calico_sample_id != "No_calico")

clone_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/"      
clone_df = read_tsv(str_interp("${clone_dir}/Consensus_genetic_clones_2023-10-05.tsv"))    

output_dir = str_interp("${output_dir}/${case_id}")
if (!dir.exists(output_dir)) {dir.create(output_dir)}

###################################################
################# Read InferCNV ###################
###################################################

#### Test: modify infercnv default plotting
for (sample_id in sample_id_list) {
    print(sample_id)
    infercnv_dir = str_interp("/diskmnt/Projects/Users/congma/outputs/infercnv_test/${setting}/${sample_id}/")
    infercnv_obj_path = str_interp("${infercnv_dir}/run.final.infercnv_obj")
    if (file.exists(str_interp("${output_dir}/${sample_id}_infercnv_ungroup.observations.txt"))) {next}

    if(!file.exists(infercnv_obj_path)) {
	    infercnv_obj_path = str_interp("${infercnv_dir}/preliminary.infercnv_obj")
    }
    infercnv_obj = readRDS(infercnv_obj_path)
    reorder_obj = infercnv_obj
    reorder_obj@tumor_subclusters = NULL
    plot_cnv(reorder_obj, 
	     out_dir = output_dir, 
	     output_filename = str_interp("${sample_id}_infercnv_ungroup"),
	     cluster_by_groups = F,
	     write_expr_matrix = T,
             x.center = mean(reorder_obj@expr.data),
	     )
}

##################################################
############ Redo InferCNV plotting ##############
##################################################

case_infercnv_df = data.frame()
case_subclone_df = data.frame()
for (sample_id in sample_id_list) {
    print(sample_id)
    infercnv_mat = data.table::fread(str_interp("${output_dir}/${sample_id}_infercnv_ungroup.observations.txt")) 
    infercnv_df = infercnv_mat %>% column_to_rownames("V1") %>% t %>% as.data.frame %>% 
	    rownames_to_column("barcode") %>% mutate(barcode=gsub('\\.',"-",barcode)) %>% 
	    mutate(sample_id=sample_id)
    case_infercnv_df = bind_rows(case_infercnv_df, infercnv_df)

    subclone_tbl = read_subclone_tbl(sample_id) %>% 
	    left_join(clone_df %>% filter(sample_id==!!sample_id) %>% mutate_all(as.character) %>% select(-sample_id)) 
    
    case_subclone_df = bind_rows(case_subclone_df, 
				 subclone_tbl %>% mutate(sample_id=sample_id))
}

gene_order_file = read.delim(file = "/diskmnt/Projects/Users/cliu/metnet_mCRC/scripts/InferCNV/gencode_v32_gene_name.txt", header=FALSE, stringsAsFactors = FALSE, sep="\t")

case_infercnv_mat = case_infercnv_df %>%
	unite(sample_barcode, c("sample_id","barcode")) %>% 
	column_to_rownames("sample_barcode") %>% 
	as.matrix
gene_order_df = gene_order_file %>% 
	filter(V1 %in% colnames(case_infercnv_mat), grepl("^chr",V2)) 
gene_order = gene_order_df$V1
chr = gene_order_df$V2
chr_level = gsub("chr","",chr) %>% as.numeric %>% unique %>% sort %>% paste0("chr",.)
chr = factor(chr, levels=chr_level)

case_infercnv_mat[is.na(case_infercnv_mat)] = 1
case_infercnv_mat = case_infercnv_mat[,gene_order]
case_infercnv_mat %>% dim

row_anno_df = case_subclone_df %>%
	mutate(sample_barcode = paste0(sample_id,"_",barcode)) %>% 
	column_to_rownames("sample_barcode") %>%
	mutate(clone_id= LETTERS[as.numeric(gsub("clone_","",genetic_clone))],
	       piece_id = sample_id) %>%
	select(clone_id, piece_id) %>% 
	.[rownames(case_infercnv_mat),]

piece_col = rainbow_hcl(length(unique(row_anno_df$piece_id))) %>% `names<-`(sort(unique(row_anno_df$piece_id)))
clone_col = rainbow_hcl(length(unique(row_anno_df$clone_id))) %>% `names<-`(sort(unique(row_anno_df$clone_id)))

max_abs_cnv = max(abs(case_infercnv_mat-1))
ht_infercnv = Heatmap(case_infercnv_mat, name = "CNV", 
        use_raster=TRUE, raster_quality=0.5,raster_resize_mat = T,
	col = colorRamp2(c(1-max_abs_cnv, 1, 1+max_abs_cnv), hcl_palette = "RdBu",reverse=T),
        column_split = chr, 
	cluster_columns = F, cluster_column_slices=F,, # show_row_dend = TRUE,
        row_split = row_anno_df$clone_id, 
        cluster_row_slices = F, 
        show_row_names = F, show_column_names = F,
        row_gap = unit(0.1, "points"),
        column_title_gp = gpar(fontsize = 10), border = TRUE,
        column_gap = unit(0, "points"),
        column_title = ifelse(1:length(chr_level) %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
        row_title = "Spatial\n(InferCNV)",
        row_title_gp = gpar(fontsize = 8), row_title_rot = 0,
        height = unit(4, "cm"),
        left_annotation = rowAnnotation(df=row_anno_df, show_annotation_name = FALSE,
            col = list("piece_id"=piece_col, "clone_id"=clone_col) 
            ),
        heatmap_legend_param = list(direction = "vertical", title_position = "topleft")
        )   

pdf(str_interp("${output_dir}/Heatmap_${case_id}_InferCNV.pdf"), w=12, h=6)
draw(ht_infercnv, merge_legend=T)
dev.off()

