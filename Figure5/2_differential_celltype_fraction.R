# RCTD subclone TME composition summary
# Jingxian Clara Liu
# 2023/01/20

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

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")
selected_tbl = read_tracking_sheet()

mode = "Sample_ref"
output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/13-RCTD/13_3-cohort/${mode}/")
celltype_order = c("Tumor", cell_types_myeloid, cell_types_lymphoid, cell_types_stroma, cell_types_breast, cell_types_liver, cell_types_pancreas)

### Set up sample to be processed -----------------------------------
region_size_df = read_tsv(REGION_SUMMARY_PTH) 

clone_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/"
clone_df = read_tsv(str_interp("${clone_dir}/Consensus_genetic_clones_2023-10-05.tsv")) %>% 
	mutate(clone_id = LETTERS[as.numeric(gsub("clone_","",genetic_clone))])
sample_id_list = clone_df$sample_id %>% unique 

sample_id_list = selected_tbl %>% filter(Study_cohort=="Discovery") %>% pull(Sample) %>% 
	intersect(sample_id_list)

### Read all RCTD composition --------------------------------------------
all_rctd_df = read_tsv(str_interp("${output_dir}/RCTD_all_celltype_cohort_2023-10-05.tsv"))
case_id_by_cancertype =  selected_tbl %>%
	filter(Sample %in% all_rctd_df$sample_id) %>% 
	arrange(cancer_type) %>% 
	pull(case_id) %>% unique

df_raw = all_rctd_df %>%
    filter(sample_id!="HT308B1-S2H5Fc2U1Z1Bs1") %>% 
    mutate_at(vars(starts_with("RCTD_")), function(x){return(ifelse(is.na(x),0,x))}) %>% 
#    select(-boundary_type) %>% 
    pivot_longer(cols=starts_with("RCTD_"), names_to="celltype", values_to="RCTD_weight") %>%
    mutate(celltype = gsub("RCTD_","",celltype)) %>%
    mutate(celltype = factor(celltype, levels=celltype_order)) %>%
    filter(!is.na(RCTD_weight))

all_celltypes = unique(df_raw$celltype)
immune_celltypes = c(cell_types_myeloid, cell_types_lymphoid)
normal_epi_celltypes = c(cell_types_liver, cell_types_pancreas, cell_types_breast)

df_toplot = df_raw %>% 
    filter(RCTD_weight>0) %>% 
    left_join(clone_df, by=c("sample_id","Filtered_tumor_regions")) %>%
    left_join(region_size_df, by=c("sample_id","Filtered_tumor_regions")) %>%
    mutate(tumor_TME=ifelse(is.na(manual_review)&Filtered_tumor_regions!=0, "Tumor", "TME")) %>% 
    mutate(celltype_class = case_when(
	celltype == "Tumor" ~ "Tumor",
	celltype %in% immune_celltypes ~ "Immune",
	celltype %in% normal_epi_celltypes ~ "Normal Epi",
	TRUE ~ "Stroma"
	)) %>% 
    mutate(celltype_class = factor(celltype_class, levels=names(celltype_class.color)))
stopifnot(nrow(df_raw %>% filter(RCTD_weight>0)) == nrow(df_toplot))

table(df_toplot %>% select(tumor_TME, celltype))

# -------------- Clone: barplot -------------
boundary = ""
boundary = "_boundary"
for (setting in c("All", "NonTumor","Immune","Lymphoid","Stroma")) {
    df_clone_toplot = df_toplot %>% filter(!is.na(clone_id))
    if (boundary=="_boundary") {
	# df_clone_toplot = df_clone_toplot %>% filter(!is.na(boundary_type))
	df_clone_toplot = df_clone_toplot %>% filter(boundary_type=="Tumor boundary")
    }
#### Fraction out of all non-tumor 
#    if (setting != "All") {df_clone_toplot = df_clone_toplot %>% filter(celltype!="Tumor")}
    df_clone_toplot = df_clone_toplot  %>%
	group_by(clone_id, case_id, celltype) %>% 
	summarise(RCTD_weight=sum(RCTD_weight, na.rm=T)) %>% 
	ungroup %>% group_by(clone_id, case_id) %>%
	mutate(TME_RCTD_sum=sum(RCTD_weight, na.rm=T)) %>%
	ungroup %>% 
	mutate(RCTD_weight = ifelse(TME_RCTD_sum==0,0,RCTD_weight/TME_RCTD_sum)) %>% 
	mutate(case_id = factor(case_id, levels=case_id_by_cancertype))
    if (setting == "NonTumor") {
	df_clone_toplot = df_clone_toplot %>% filter(celltype!="Tumor")
    } else if (setting == "Immune") {
	df_clone_toplot = df_clone_toplot %>% 
		filter(celltype %in% c(cell_types_lymphoid, cell_types_myeloid))
    } else if (setting == "Lymphoid") {
	df_clone_toplot = df_clone_toplot %>% 
		filter(celltype %in% c(cell_types_lymphoid))
    } else if (setting == "Stroma") {
	df_clone_toplot = df_clone_toplot %>% 
		filter(celltype %in% c(cell_types_stroma, cell_types_breast, cell_types_liver, cell_types_pancreas))
    }
    plot_title = str_interp("${setting} Cell Type Fraction")
    if (boundary=="_boundary") {plot_title=str_interp("${plot_title} (boundary only)")}
    p = ggplot(df_clone_toplot, aes(x=clone_id, y=RCTD_weight*100, fill=celltype)) +
	geom_bar(stat="identity") +
	facet_grid(~case_id, space="free", scales="free_x") +
	scale_fill_manual(values=celltype.color[unique(as.character(df_clone_toplot$celltype))]) +
	scale_y_continuous(guide = guide_axis(n.dodge = 2), breaks=scales::pretty_breaks(),
			   expand=expansion(c(0,0.1))) + 
	theme_bw_nogrid +
	theme(legend.position="bottom", axis.text=element_text(color="black", size=12)) +
	theme(legend.text = element_text(size=3), legend.key.size = unit(2, 'mm')) +
	guides(fill = guide_legend(nrow = 2)) +
	labs(y="Cell type percentage", x="Case / Clone", fill="Cell Type", 
	     title=plot_title)
    ggsave(str_interp("${output_dir}/BySubclone/Barplot_Tumor_subclones_RCTD_weights_${setting}${boundary}.pdf"), plot=p, w=7, h=3.5)
}


# -------------- Clone x cell type: test -------------
barcode_rctd_sum = df_toplot %>% group_by(sample_id,barcode) %>% summarise(sum=sum(RCTD_weight)) %>% pull(sum)
stopifnot(all.equal(barcode_rctd_sum,rep(1, length(barcode_rctd_sum))))
test_df = data.frame()
for (case_id in unique(df_toplot$case_id)) {
    for (celltype in unique(df_toplot$celltype)) {
	df_test = df_toplot %>% 
	    filter(celltype==!!celltype, case_id==!!case_id) %>% 
	    filter(clone_id != "TME")
        if (nrow(df_test)==0) {next}
	if (length(unique(df_test$clone_id))==1) {next}
	test = pairwise.wilcox.test(df_test$RCTD_weight, df_test$clone_id, p.adjust.method='fdr')
	test_df = test_df %>% 
	    bind_rows(tidy(test) %>% mutate(case_id=case_id, celltype=celltype))
    }
}
test_df = test_df %>% arrange(p.value)
# test_df %>% write_tsv(str_interp("${output_dir}/BySubclone/wilcoxon_test_rctd_subclone.tsv"))

##### scatter plot: 
test_df = read_tsv(str_interp("${output_dir}/BySubclone/wilcoxon_test_rctd_subclone.tsv"))
test_df_toplot = test_df %>%
    filter(celltype!="Tumor") %>%
    mutate(FDR=p.adjust(p.value)) %>% 
    left_join(selected_tbl %>% distinct(case_id, cancer_type)) %>% 
    group_by(cancer_type) %>%
    arrange(FDR, p.value) %>%
    mutate(rank=row_number()) %>%
    ungroup %>% 
#    mutate(label=ifelse(p.value<0.05 & rank<=5, case_id, NA_character_)) %>% 
    mutate(label=paste0(case_id," (",celltype,")")) %>%
#    mutate(label=case_id) %>% 
    filter(p.value<0.05)
p = ggplot(test_df_toplot,
	   aes(x=rank, y=-log10(FDR), fill=celltype)) +
	facet_grid(. ~ cancer_type, scales="free", space="free") +
	geom_bar(stat="identity") +
#	geom_text_repel(aes(label=label), size=2, max.overlaps=Inf, min.segment.length = 0) +
	geom_text(aes(label=label, color=celltype), size=2, angle=90, hjust=-0.1) +
#	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	scale_fill_manual(values=celltype.color) + scale_color_manual(values=celltype.color) +
	scale_y_continuous(breaks=scales::pretty_breaks(), expand=expansion(c(0,0.25))) +
	scale_x_continuous(expand=expansion(c(0,0))) +
	theme_bw_nogrid +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	theme(legend.position="none", strip.background = element_rect(color=NA)) +
	labs(x="Rank", title="Differential TME infiltration between spatial subclones")
ggsave(str_interp("${output_dir}/BySubclone/wilcoxon_test_rctd_subclone_scatter.pdf"), 
       plot=p, w=6, h=2.5)

##### Example: HT268B1
case_id = "HT268B1"
celltype="Macrophage"

case_id = "HT397B1"
celltype="Macrophage"
celltype="T"

case_id = "HT308B1"
celltype="Macrophage"

example_df = df_toplot %>% 
    filter(case_id==!!case_id, celltype==!!celltype) %>% 
    filter(clone_id!="TME")
all_clones = unique(example_df$clone_id)
all_clones_combo = combn(all_clones, 2) %>% as.data.frame %>% as.list 
p = ggplot(example_df, aes(x=clone_id, y=RCTD_weight, color=celltype)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.2, size=0.1) +
    stat_compare_means(comparisons=all_clones_combo, size=3, aes(label="..p.format..")) +
    scale_color_manual(values=celltype.color) +
    scale_y_continuous(expand=expansion(c(0,0.1)), breaks=scales::pretty_breaks()) +
    theme_bw_nogrid +
    theme(legend.position="none", axis.text=element_text(color="black"),
	  plot.title = element_text(size=10)) +
    labs(title=str_interp("${case_id}\n(${celltype})"), x="Clone ID", y="Cell Type Fraction", color="Cell Type")
ggsave(str_interp("${output_dir}/BySubclone/wilcoxon_test_rctd_subclone_example_${case_id}_${celltype}.pdf"), 
       plot=p, w=1.5, h=2.5)

##### Example
