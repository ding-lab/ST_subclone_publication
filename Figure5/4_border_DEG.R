# Calculate genes differentiall expressed in the tumor boundary regions
# Jingxian Clara Liu
# 2023/01/19

library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(optparse)
library(googlesheets4)
library(ComplexHeatmap)
# library(ggpubr)

# Test parameters
ASSAY = "Spatial"
SLOT = "data"
# ASSAY = "SCT"
# SLOT = "scale.data"
WIDTH = 2
TOP_N=10
MIN_SPOT=10

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")

sample_id = "HT112C1-U1_ST_Bn1"
interaction_path = "/diskmnt/Projects/Users/cliu/metnet_mCRC/scripts/cellphonedb/cpdb_interaction_input_gene_name.tsv"
output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/16-CCI/16_3-boundary_DEG/width_${WIDTH}/")

selected_tbl = read_tracking_sheet()

### Parameters ------------------------------------------------------
option_list = list(
  make_option(c("-s", "--sample_id"),
        type="character",
        default=NULL,
        help="Sample name for ST",
        metavar="character"),
  make_option(c("-o", "--st_obj_path"),
        type="character",
        default="",
        help="Path to ST seurat object",
        metavar="character"),
  make_option(c("-a", "--subclone_path"),
        type="character",
        default="",
        help="Path to subclone annotation file",
        metavar="character"),
  make_option(c("-d", "--output_dir"),
        type="character",
        default=output_dir,
        help="Directory to output",
        metavar="character"),
  make_option(c("-r", "--interaction_path"),
        type="character",
        default=interaction_path,
        help="Path to gene-gene interaction file",
        metavar="character"),
  make_option(c("-t", "--top_n"),
        type="integer",
        default=10,
        help="Number of top interactions to plot",
        metavar="integer"),
  make_option(c("-w", "--width"),
        type="integer",
        default=4,
        help="Width of boundary region in number of spots",
        metavar="integer")
  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

sample_id=opt$sample_id
st_obj_path=opt$st_obj_path
subclone_path=opt$subclone_path
output_dir=opt$output_dir
interaction_path=opt$interaction_path
TOP_N=opt$top_n
WIDTH=opt$width

### Set up sample to be processed -----------------------------------
output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/16-CCI/16_3-boundary_DEG/width_${WIDTH}/")
if(!dir.exists(output_dir)) {dir.create(output_dir)}
sample_output_dir = str_interp("${output_dir}/${sample_id}/")
if(dir.exists(sample_output_dir)) {unlink(sample_output_dir, recursive=T)}
if(!dir.exists(sample_output_dir)) {dir.create(sample_output_dir)}

subclone_tbl_toadd = read_subclone_tbl(sample_id)
st = read_st_obj(sample_id)

st = NormalizeData(st, verbose=F, assay="Spatial")
st = ScaleData(st, assay="Spatial")
st <- AddMetaData(obj=st, metadata=subclone_tbl_toadd)

interaction_tbl = read_tsv(interaction_path) %>% filter(gene_name_a!=gene_name_b)

if (WIDTH==1) {
  boundary_region_cells = st$TME_boundary %>% .[!is.na(.)] %>% names
} else if (WIDTH==2) {
  boundary_region_cells = c(
    st$Tumor_boundary %>% .[!is.na(.)] %>% names,
    st$TME_boundary %>% .[!is.na(.)] %>% names
  )
} else if (WIDTH==4) {
  boundary_region_cells = c(
    st$Tumor_boundary %>% .[!is.na(.)] %>% names,
    st$TME_boundary %>% .[!is.na(.)] %>% names,
    st$Tumor_before_boundary %>% .[!is.na(.)] %>% names,
    st$TME_before_boundary %>% .[!is.na(.)] %>% names
  )
}
  
width4_boundary_region_cells = c(
    st$Tumor_boundary %>% .[!is.na(.)] %>% names,
    st$TME_boundary %>% .[!is.na(.)] %>% names,
    st$Tumor_before_boundary %>% .[!is.na(.)] %>% names,
    st$TME_before_boundary %>% .[!is.na(.)] %>% names
  )

if (SLOT=="scale.data") { 
  st_exp_mat = st[[ASSAY]]@scale.data
}
if (SLOT=="data") {
  st_exp_mat= st[[ASSAY]]@data 
}
DefaultAssay(st) = ASSAY

# Yes boundary / No boundary
st$region_type = case_when(
  Cells(st) %in% boundary_region_cells ~ "Boundary",
  st$Filtered_tumor_regions != 0 ~ "Tumor",
  st$Filtered_tumor_regions == 0 ~ "TME",
  TRUE ~ NA_character_ 
) %>% factor(., levels=c("Tumor","Boundary","TME"))

Idents(st) = 'region_type'
bound_markers = FindMarkers(st, ident.1="Boundary")
bound_tumor_markers = FindMarkers(st, ident.1="Boundary", ident.2="Tumor")
bound_tme_markers = FindMarkers(st, ident.1="Boundary", ident.2="TME")

markers = bind_rows(
  bound_markers %>% mutate(cluster="Boundary") %>% rownames_to_column("gene"),
  bound_tumor_markers %>% mutate(cluster="Boundary_vs_Tumor") %>% rownames_to_column("gene"),
  bound_tme_markers %>% mutate(cluster="Boundary_vs_TME") %>% rownames_to_column("gene")
  )

# Yes boundary / leave out middle boundary / No boundary
st$region_type_leaveout = case_when(
  Cells(st) %in% boundary_region_cells ~ "Boundary",
  Cells(st) %in% width4_boundary_region_cells ~ "Leaveout",
  st$Filtered_tumor_regions != 0 ~ "Tumor",
  st$Filtered_tumor_regions == 0 ~ "TME",
  TRUE ~ NA_character_ 
) 

st_lo = subset(st, region_type_leaveout != "Leaveout") 
st_lo$region_type_leaveout = factor(st_lo$region_type_leaveout, 
				    levels=c("Tumor","Boundary","TME"))

Idents(st_lo) = 'region_type_leaveout'
bound_markers_lo = FindMarkers(st_lo, ident.1="Boundary")
bound_tumor_markers_lo = FindMarkers(st_lo, ident.1="Boundary", ident.2="Tumor")
bound_tme_markers_lo = FindMarkers(st_lo, ident.1="Boundary", ident.2="TME")

markers_lo = bind_rows(
  bound_markers_lo %>% mutate(cluster="Boundary_lo") %>% rownames_to_column("gene"),
  bound_tumor_markers_lo %>% mutate(cluster="Boundary_vs_Tumor_lo") %>% rownames_to_column("gene"),
  bound_tme_markers_lo %>% mutate(cluster="Boundary_vs_TME_lo") %>% rownames_to_column("gene")
  )

bind_rows(markers, markers_lo) %>%  
	write_tsv(str_interp("${sample_output_dir}/${sample_id}_boundary_DEG.tsv"))  

plot_mode="lo"
for (plot_mode in c("all","lo")) {
selected_markers = bind_rows(markers, markers_lo) %>%
	mutate(mode = ifelse(grepl("_lo$",cluster), "lo", "all")) %>%
	filter(mode==plot_mode) %>%
	mutate(cluster = gsub("_lo$", "", cluster)) %>%
	pivot_wider(id_cols="gene", names_from="cluster", values_from="avg_log2FC") %>%
	filter_at(vars(starts_with("Boundary")), all_vars(!is.na(.)&.>0)) %>%
	arrange(desc(Boundary)) %>% slice_head(n=20) %>% pull(gene)

markers_toplot_df = markers %>%
#	filter(pct.1-pct.2 > 0.2) %>%
	filter(avg_log2FC > 0) %>% 
	arrange(desc(avg_log2FC), pct.1-pct.2) %>% 
	group_by(cluster) %>%
	slice_head(n=20) 
markers_toplot = markers_toplot_df$gene %>% unique
markers_toplot = c(selected_markers, setdiff(markers_toplot, selected_markers))

mat = st_exp_mat[markers_toplot,] %>% as.matrix %>% t %>% scale %>% t
pos_95_mat = quantile(mat[which(mat>0)], 0.95, na.rm=T)
neg_95_mat = quantile(mat[which(mat<0)], 0.05, na.rm=T)

if (plot_mode=="all") {
  col_split = st$region_type
  df_anno = data.frame(`region_type`=st$region_type)
  ha_col = HeatmapAnnotation(
	df = df_anno, 
	col = list(`region_type`=t_b_tme.color)
  )
} else if (plot_mode=="lo") {
  col_split = st$region_type_leaveout %>% factor(., levels=c("Tumor","Boundary","TME","Leaveout"))
  df_anno = data.frame(`region_type`=st$region_type_leaveout)
  ha_col = HeatmapAnnotation(
	df = df_anno,
	col = list(`region_type`=c(t_b_tme.color, c("Leaveout"=NA)))
  )
}

ht = Heatmap(
	mat,
	name = "Normalized\nExpression\nZ score",
	col = circlize::colorRamp2(c(neg_95_mat, 0, pos_95_mat), hcl_palette = "RdBu", reverse = T),
	show_column_names=F,
	column_split = col_split, 
	row_split = ifelse(markers_toplot%in%selected_markers, "Both", ""),
	cluster_row_slices=F, cluster_column_slices=F,
	cluster_rows=T, cluster_columns=T,
	row_names_gp = gpar(fontsize = 8),

	show_row_dend = FALSE,
	show_column_dend = FALSE,

	row_names_side = "left",
	)

pdf(str_interp("${sample_output_dir}/${sample_id}_boundary_DEG_${plot_mode}.pdf"), w=12, h=7)
draw(ha_col %v% ht)
dev.off()

df_toplot = st_exp_mat[markers_toplot,] %>% 
	as.data.frame %>%
	rownames_to_column("gene") %>%
	pivot_longer(
	    cols = !gene, 
	    names_to = "barcode",
	    values_to = "expression"
	) %>% 
	left_join(df_anno %>% rownames_to_column("barcode")) %>% 
	filter(region_type!="Leaveout", gene %in% selected_markers) %>% 
	mutate(region_type = factor(region_type, levels=c("Tumor","Boundary","TME")))

p = ggplot(df_toplot, aes(x=region_type, y=expression)) +
	geom_violin(aes(fill=region_type)) +
	facet_wrap(~ gene, scale="free_y") +
	theme_bw() +
	scale_fill_manual(values=t_b_tme.color) +
	labs(x="Region type", y="Expression", fill="Region type")

ggsave(str_interp("${sample_output_dir}/${sample_id}_boundary_DEG_${plot_mode}_violin.pdf"), 
       plot=p, w=10, h=5)

}
