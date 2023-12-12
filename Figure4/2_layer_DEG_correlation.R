# Layer-gene correation test to identify gene expression associated with tumor depth
# Jingxian Clara Liu
# 2023/10/03

library(Matrix, lib.loc="/diskmnt/Projects/Users/cliu/software/anaconda3/envs/r41-base/lib/R/library")
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(optparse)
library(umap)
library(googlesheets4)
library(ComplexHeatmap)
library(colorspace)
library(ggpubr)
library(rstatix)
library(patchwork)
library(ggpmisc)
library(future)
library(circlize)
library(clusterProfiler, lib.loc="/diskmnt/Projects/Users/cliu/software/anaconda3/envs/r4-base/lib/R/library")
library(msigdbr, lib.loc="/diskmnt/Projects/Users/cliu/software/anaconda3/envs/r4-base/lib/R/library")

plan("multicore", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2) # 4GB

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")
selected_tbl = read_tracking_sheet()

sample_id = "HT260C1-Th1K1Fc2U1Z1Bs1"
sample_id = "HT268B1-Th1K3Fc2U1Z1Bs1"
sample_id = "HT270P1-S1H1Fs5U1Bp1"
sample_id = "HT206B1-S1Fc1U2Z1B1"

FC = "FC5"
# comp_type = "Layer_num"
comp_type = "Layer_frac"

covar = "none"
covar = "purity"

### Parameters ------------------------------------------------------
option_list = list(
  make_option(c("-s", "--sample_id"),
	type="character",
	default=NULL,
	help="Sample name for ST",
	metavar="character"),
  make_option(c("-c", "--covar"),
	type="character",
	default="none",
	help="Covariate to include",
	metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

sample_id=opt$sample_id
covar=opt$covar

### Read ST obj -------------------------------------------------
if (covar == "none") {
    output_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/40-Layer_DEG/40_1-Layer_DEG/"
} else if (covar == "purity") {
    output_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/40-Layer_DEG/40_3-Layer_DEG_purity/"
    comp_type = "Layer_frac"
}

sample_output_dir = str_interp("${output_dir}/${sample_id}")
if (!dir.exists(sample_output_dir)) {dir.create(sample_output_dir)}
    
subclone_tbl = read_subclone_tbl(sample_id) %>%
    mutate(FC5 = ifelse(FC5=="unknown", NA_real_,as.numeric(FC5))) %>% 
    mutate_at(vars(starts_with("FC")), function(x){ifelse(x==0, NA_real_,x)})
st = read_st_obj(sample_id)
st = AddMetaData(st, subclone_tbl[Cells(st),])

region_big = subclone_tbl %>% 
    filter(Filtered_tumor_regions!="0") %>%
    group_by(Filtered_tumor_regions) %>% 
    summarise(n_spot=n(),
	      FC1_max=max(FC1, na.rm=T), FC3_max=max(FC3, na.rm=T),
	      FC5_max=max(FC5, na.rm=T)) %>% 
    filter(n_spot>=50, (FC1_max>=3|FC3_max>=3|FC5_max>=3)) %>% 
    pull(Filtered_tumor_regions) %>% as.character
stopifnot(length(region_big)>0)

st_tumor = subset(st, Filtered_tumor_regions != "0")
st_tumor_big = subset(st, Filtered_tumor_regions %in% region_big)
st_tumor_big$FC = st_tumor_big@meta.data[,FC]
st_tumor_big$FC_max = st_tumor_big@meta.data %>% group_by(Filtered_tumor_regions) %>% 
    mutate(FC_max = max(FC, na.rm=T)) %>% pull(FC_max)
st_tumor_big$FC_frac = (st_tumor_big$FC-1) / (st_tumor_big$FC_max-1)

DefaultAssay(st_tumor_big) = "Spatial"
st_tumor_big = FindVariableFeatures(st_tumor_big)
st_tumor_big = NormalizeData(st_tumor_big)
st_tumor_big = ScaleData(st_tumor_big)
st_exp_mat = st_tumor_big[["Spatial"]]@data

### Plot results  -------------------------------------------------
cor_adj_df = read_tsv(str_interp("${sample_output_dir}/${sample_id}_layer_${FC}_DEG.tsv")) %>% 
	filter(comp_type==!!comp_type, !is.na(rho))
if (covar == "none") {
    tumor_normal_deg = read_tsv(str_interp("${sample_output_dir}/${sample_id}_TumorNormal_DEG.tsv"))
    cor_adj_df = cor_adj_df %>% left_join(tumor_normal_deg, by="gene") 
}
cor_adj_df = cor_adj_df %>% arrange(desc(abs(rho)))

cor_adj_df_dir = cor_adj_df %>% 
	mutate(direction = ifelse(rho>0, "Center", "Peri")) %>% 
	mutate(gene_label=ifelse(p.value<1e-5, gene, ""))
# cor_adj_df_dir %>% filter(avg_log2FC>0, rho<0)  

####### 1. Selected genes: cor & spatial plot --------
center_deg = cor_adj_df %>% 
	filter(rho>0) %>% 
	filter(!grepl("^MT|^RP", gene))
top_layer_deg = center_deg %>% pull(gene) %>% unique %>% head(n=3)
periphery_deg = cor_adj_df %>% 
	filter(rho<0) %>% 
	filter(!grepl("^MT|^RP", gene)) 
bot_layer_deg = periphery_deg %>% pull(gene) %>% unique %>% head(n=3) 
deg_toplot = c(top_layer_deg, bot_layer_deg)
# if (grepl("HT260", sample_id)) {top_layer_deg = c("VEGFA", "TFF3", "DDIT4")}
if (sample_id=="HT206B1-S1Fc1U2Z1B1") {deg_toplot=c("NDRG1","TGFBI","CA9","TUBA1B","NDUFA4","TOMM40")}
deg_toplot

x_var = ifelse(comp_type == "Layer_frac", "FC_frac", "FC")
st_exp_df = st_exp_mat %>% as.data.frame %>% 
        rownames_to_column("gene") %>%
        pivot_longer(cols=!gene, names_to="barcode", values_to="expr") %>% 
        left_join(st_tumor_big@meta.data %>% select(barcode, Filtered_tumor_regions, !!x_var)) %>%
        filter(gene %in% deg_toplot) %>% 
	mutate(gene_type = ifelse(gene %in% center_deg$gene, "Center", "Periphery")) %>%
        mutate(gene=factor(gene,levels=unique(deg_toplot))) 
   
p = ggplot(st_exp_df, aes_string(x=x_var, y="expr")) +
	geom_boxplot(aes_string(group=x_var), outlier.shape=NA) +
	# geom_point(size=0.1) +
        facet_wrap(~gene, nrow=1, scale='free_y') +
        theme_bw() +
        stat_cor(aes(color=gene_type), p.accuracy=0.001) +
        stat_poly_line(aes(color=gene_type)) +
	scale_color_manual(values=cntr_peri.color) +
	labs(x="Spot depth", y="Expression") +
	theme_bw_nogrid +
	scale_x_continuous(breaks= scales::pretty_breaks()) +
	theme(strip.background = element_rect(fill=NA, color=NA),
	      axis.text = element_text(size=12),
	      strip.text = element_text(face="bold", size=12), 
	      legend.position="none")
p_st = SpatialFeaturePlot(st_tumor, deg_toplot, stroke=0, crop=F, ncol=length(deg_toplot)) &
	theme(legend.position="none")
ggsave(str_interp("${sample_output_dir}/${sample_id}_${comp_type}_DEG_most_sig_selected.pdf"), 
       plot=wrap_plots(p, p_st, ncol=1, heights=c(1,1.5)), w=2*length(deg_toplot), h=5)

####### 2. sample pathway -------
# Hallmark
H_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
# C5: GOBP
C5BP_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5") %>% 
    dplyr::select(gs_name, human_gene_symbol) %>% 
    filter(grepl("^GOBP_",gs_name))

#### Tumor region gene (avg_log2FC>0), higher in center (rho>0), sig layer (p.val.adj<0.05)
tumor_gene_cent = cor_adj_df %>%
    filter(rho>0.1, p.val.adj<0.05) %>%
#    slice_max(abs(rho), n=150) %>%
    pull(gene) %>% unique
length(tumor_gene_cent)
tumor_gene_cent_path = enricher(tumor_gene_cent, TERM2GENE=bind_rows(H_t2g, C5BP_t2g))
tumor_gene_cent_path %>% as.data.frame %>% write_tsv(str_interp("${sample_output_dir}/${comp_type}_${FC}_tumor_high_center_pathway_H_C5BP.tsv"))
# tumor_gene_cent_path %>% as.data.frame %>% filter(grepl("HYPOXIA", Description))

tumor_gene_peri = cor_adj_df %>%
    filter(rho<(-0.1), p.val.adj<0.05) %>% 
#    slice_max(abs(rho), n=150) %>%
    pull(gene) %>% unique
length(tumor_gene_peri)
tumor_gene_peri_path = enricher(tumor_gene_peri, TERM2GENE=bind_rows(H_t2g, C5BP_t2g))
tumor_gene_peri_path %>% as.data.frame %>% write_tsv(str_interp("${sample_output_dir}/${comp_type}_${FC}_tumor_high_periphery_pathway_H_C5BP.tsv"))

#### TME region gene (avg_log2FC<0), higher in center (rho>0), sig layer (p.val.adj<0.05)
if (F) {
tme_gene_cent = cor_adj_df %>%
    filter(avg_log2FC<0, rho>0.1, p.val.adj<0.05) %>% 
    pull(gene) %>% unique
length(tme_gene_cent)
tme_gene_cent_path = enricher(tme_gene_cent, TERM2GENE=bind_rows(H_t2g, C5BP_t2g))
tme_gene_cent_path %>% as.data.frame %>% write_tsv(str_interp("${sample_output_dir}/${comp_type}_${FC}_tme_high_center_pathway_H_C5BP.tsv"))

tme_gene_peri = cor_adj_df %>%
    filter(avg_log2FC<0, rho<(-0.1), p.val.adj<0.05) %>% 
    pull(gene) %>% unique
length(tme_gene_peri)
tme_gene_peri_path = enricher(tme_gene_peri, TERM2GENE=bind_rows(H_t2g, C5BP_t2g))
tme_gene_peri_path %>% as.data.frame %>% write_tsv(str_interp("${sample_output_dir}/${comp_type}_${FC}_tme_high_periphery_pathway_H_C5BP.tsv"))
}

all_path_df = data.frame()
for (tissue in c("tumor","tme")) {
    for (layer_dir in c("center","periphery")) {
	dir_color = cntr_peri.color[str_to_title(layer_dir)] 
	path_df_path = str_interp("${sample_output_dir}/${comp_type}_${FC}_${tissue}_high_${layer_dir}_pathway_H_C5BP.tsv")
	if (!file.exists(path_df_path)) {next}
	path_df = read_tsv(path_df_path)
	all_path_df = bind_rows(all_path_df, path_df %>% mutate(dir=layer_dir))
	if (nrow(path_df) == 0) {next}
        path_toplot = path_df %>%
		filter(grepl("HALLMARK", ID)) %>%
		mutate(ID = gsub("_","_",gsub("HALLMARK_", "", ID))) %>% 
		mutate(p.adjust = p.adjust(pvalue,method="fdr")) %>% 
		arrange(p.adjust) %>% slice_head(n=20)
	if (nrow(path_toplot) == 0) {next}
	path_toplot$ID = factor(path_toplot$ID, levels=rev(path_toplot$ID))
	p = ggplot(path_toplot, aes(x=-log10(p.adjust), y=ID)) +
		geom_bar(stat="identity", fill=dir_color) +
		scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
		geom_vline(xintercept=-log10(0.05), linetype="dashed", size=1) +
		theme_bw_nogrid +
		theme(axis.text = element_text(size=12, color="black"), 
		      axis.title = element_text(size=12, color="black"), 
		      legend.position="bottom") +
		labs(title=str_interp("Enriched in ${layer_dir}"), y="", x="-log10(FDR)")
	ggsave(str_interp("${sample_output_dir}/${comp_type}_${FC}_${tissue}_high_${layer_dir}_pathway_H_C5BP.pdf"), plot=p, w=5, h=4)
    }
}
if (nrow(all_path_df)>0) {
    path_toplot = all_path_df %>% 
	mutate(ID = gsub("_"," ",gsub("HALLMARK_|GOBP_", "", ID))) %>% 
	mutate(ID = str_wrap(ID, width=40)) %>%
	group_by(dir) %>% 
	slice_min(p.adjust, n=15) %>% 
	mutate(rank=ifelse(dir=="center", -log10(p.adjust), log10(p.adjust))) %>%
	mutate(dir = str_to_title(dir)) %>%
	arrange(desc(rank))
    path_toplot$ID = factor(path_toplot$ID, levels=rev(path_toplot$ID))
    p = ggplot(path_toplot, aes(x="0", y=ID)) +
#	geom_bar(aes(fill=dir), stat="identity") +
	geom_point(aes(size=-log10(p.adjust), fill=dir), shape=21) +
#	geom_vline(xintercept=-log10(0.05), linetype="dashed", size=1) +
	scale_fill_manual(values=cntr_peri.color) +
	theme_bw_nogrid +
	theme(axis.text = element_text(size=8, color="black"),
	      axis.text.x = element_blank(),
	      axis.title = element_text(size=12, color="black"),
	      legend.position="bottom") +
	scale_x_discrete(expand=c(0,0)) +
#	scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
	labs(title="", y="", x="")
    ggsave(str_interp("${sample_output_dir}/${comp_type}_${FC}_both_pathway_H_C5BP_${sample_id}.pdf"), 
	   plot=p, w=4, h=8)
}

####### 3. sn reference top genes -------
snRNA_id = read_column_val(sample_id, "Paired_snRNA_id")
snRNA_meta = read_column_val(sample_id, "Paired_snRNA_metadata") %>% 
	read_tsv %>% column_to_rownames("barcode")
sn = read_column_val(sample_id, "Paired_snRNA_obj") %>% readRDS %>% 
	AddMetaData(snRNA_meta[Cells(.),])
doublet_col = ifelse(grepl("_combo", snRNA_id), "predicted_doublet_rna", "predicted_doublet")
sn = subset(sn, cells = Cells(sn)[(!sn@meta.data[,doublet_col]) & (sn$cell_type_final!="Doublet")])
DefaultAssay(sn) = "RNA"
sn = NormalizeData(sn)
sn = FindVariableFeatures(sn)
sn = ScaleData(sn)
sn_rna = DietSeurat(sn, assays = "RNA")

top_layer_deg = cor_adj_df %>% 
	filter(rho>0) %>% 
	filter(!grepl("^MT|^RP", gene)) %>% 
	pull(gene) %>% unique %>% head(n=25)
bot_layer_deg = cor_adj_df %>% 
	filter(rho<0) %>%
#	filter(avg_log2FC>0) %>% 
	filter(!grepl("^MT|^RP", gene)) %>%
	pull(gene) %>% unique %>% head(n=25) %>% rev
top_layer_deg = tumor_gene_cent
bot_layer_deg = tumor_gene_peri
deg_toplot = c(top_layer_deg, bot_layer_deg)
deg_dir = c(rep("Center", length(top_layer_deg)), rep("Periphery", length(bot_layer_deg)))
deg_toplot

gene_celltype_mat_list = AverageExpression(sn_rna, group.by="cell_type_final", assay="RNA",
					   features = deg_toplot)
mat_toplot = gene_celltype_mat_list[["RNA"]] %>% t %>% scale %>% t

max_mat = quantile(abs(mat_toplot), 0.99, na.rm=T)
ht = Heatmap(
    mat_toplot,
    name = "AvgExpr",
    row_names_gp = gpar(fontsize = 3),
    col = circlize::colorRamp2(c(-max_mat, 0, max_mat), hcl_palette = "RdBu", reverse = T),
    cluster_rows = T, 
    cluster_columns = T,
    cluster_row_slices = F,
    row_split = deg_dir, 
    show_row_dend = F, show_column_dend = F,
    row_names_side = "left"
    )

pdf(str_interp("${sample_output_dir}/Heatmap_center_top_${FC}_${sample_id}.pdf"), w=5, h=12)
draw(ht,  heatmap_legend_side = "right", merge_legend = T, annotation_legend_side = "right")
dev.off()

