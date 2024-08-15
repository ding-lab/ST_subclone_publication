# Cohort-level pathway analysis with genes enriched in the center/periphery of tumor regions
# Jingxian Clara Liu
# 2023/10/04

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
library(clusterProfiler, lib.loc="/diskmnt/Projects/Users/cliu/software/anaconda3/envs/r4-base/lib/R/library")
library(msigdbr, lib.loc="/diskmnt/Projects/Users/cliu/software/anaconda3/envs/r4-base/lib/R/library")

plan("multicore", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2) # 4GB

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")
selected_tbl = read_tracking_sheet()

sample_id = "HT260C1-Th1K1Fc2U1Z1Bs1"
input_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/40-Layer_DEG/40_2-Layer_DEG_region/")
output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/40-Layer_DEG/40_2-Layer_DEG_region_cohort/")

### Set up sample to be processed -----------------------------------
case_id_list = selected_tbl %>% filter(Study_cohort=="Discovery") %>%
	pull(case_id) %>% sort
# case_id_list = c("HT260C1","HT397B1","HT268B1","HT270P1","HT112C1")
sample_id_list = selected_tbl %>% 
	filter(Study_cohort=="Discovery", case_id %in% case_id_list) %>%
	pull(Sample) %>% sort

all_cor_df = data.frame() 
for (sample_id in sample_id_list) {
    cor_path = str_interp("${input_dir}/${sample_id}/${sample_id}_layer_region_FC_DEG.tsv")
    if(!file.exists(cor_path)) {next}
    cor_adj_df = read_tsv(cor_path)
    all_cor_df = bind_rows(all_cor_df, cor_adj_df %>% mutate(sample_id=sample_id))
}

comp_type="FC5"

setting="any"
FDR_cutoff = 0.1
rho_abs_cutoff = 0.1

setting="all"
FDR_cutoff = 1
rho_abs_cutoff = 0

all_cor_df_select = all_cor_df %>% 
    left_join(selected_tbl %>% select(Sample, case_id, cancer_type), by=c("sample_id"="Sample")) %>% 
    filter(comp_type==!!comp_type)

all_cor_df_filter = all_cor_df_select %>% 
    mutate(up_gene_in_region = rho> rho_abs_cutoff & p.val.adj<=FDR_cutoff,
	   dn_gene_in_region = rho< (-rho_abs_cutoff) & p.val.adj<=FDR_cutoff) %>% 
    group_by(case_id, gene, cancer_type, sample_id, Filtered_tumor_regions) %>% 
    summarise(up_gene_in_sample = any(up_gene_in_region), 
	      dn_gene_in_sample = any(dn_gene_in_region), 
	      max_rho = max(rho, na.rm=T), min_rho = min(rho, na.rm=T)) %>% 
    ungroup

all_cor_df_filter %>% distinct(case_id, cancer_type) %>% arrange(cancer_type) 
# all_cor_df_filter %>% distinct(sample_id, cancer_type) %>% arrange(cancer_type)
if (setting=="all") { 
    all_cor_df_filter = all_cor_df_filter %>% 
        group_by(case_id, gene, cancer_type) %>% 
        summarise(up_gene_in_case = all(up_gene_in_sample), 
	          dn_gene_in_case = all(dn_gene_in_sample), 
	          max_rho = max(max_rho, na.rm=T), min_rho = min(min_rho, na.rm=T))
} else if (setting=="any") {
    all_cor_df_filter = all_cor_df_filter %>% 
        group_by(case_id, gene, cancer_type) %>% 
        summarise(up_gene_in_case = any(up_gene_in_sample), 
	          dn_gene_in_case = any(dn_gene_in_sample), 
	          max_rho = max(max_rho, na.rm=T), min_rho = min(min_rho, na.rm=T))
}

all_cor_df_filter = all_cor_df_filter %>% 
    group_by(gene, cancer_type) %>% 
    summarise(n_up = sum(up_gene_in_case), n_dn=sum(dn_gene_in_case),
    	      max_rho = max(max_rho, na.rm=T), min_rho = min(min_rho, na.rm=T)) %>% 
    mutate(n_up = replace_na(n_up, 0), n_dn = replace_na(n_dn, 0))
all_cor_df_filter %>% arrange(desc(n_up), desc(max_rho)) %>% as.data.frame %>% head(n=10)
all_cor_df_filter %>% arrange(desc(n_dn), min_rho) %>% as.data.frame %>% head(n=10)

for (cancer in c("all",unique(all_cor_df_filter$cancer_type))) {
    print(cancer)
    all_cor_df_plot = all_cor_df_filter
    if (cancer!="all") {
        gene_order = all_cor_df_filter %>% 
	    mutate(is_cancer = ifelse(cancer_type==!!cancer, "is_cancer", "not_cancer")) %>% 
	    group_by(gene, is_cancer) %>% 
	    summarise(n_up = sum(n_up), n_dn=sum(n_dn)) %>% 
	    pivot_wider(id_cols=c(gene), names_from="is_cancer", values_from=c(n_up, n_dn), values_fill=0) %>% 
	    mutate(n_up=n_up_is_cancer-n_up_not_cancer, 
	           n_dn=n_dn_is_cancer-n_dn_not_cancer) %>%
	    filter(n_up>=0|n_dn>=0) %>%
	    arrange(desc(n_up_is_cancer-n_dn_is_cancer)) %>% pull(gene)
#	all_cor_df_plot = all_cor_df_plot %>% filter(cancer_type==cancer)
    } else {
        gene_order =  all_cor_df_plot %>% 
	    group_by(gene) %>%
    	    summarise(n_up = sum(n_up), n_dn=sum(n_dn),
		      max_rho = max(max_rho, na.rm=T), min_rho = min(min_rho, na.rm=T)) %>%
	    arrange(desc(n_up-n_dn)) %>% pull(gene)
    }
    if (length(gene_order)>50) { gene_order =  c(head(gene_order,25), tail(gene_order,25)) }

    all_cor_df_plot = all_cor_df_plot %>%
#        filter(!grepl("^RP|^RL|^MT", gene)) %>%
#        filter(n_up>2 | n_dn>1) %>%
	filter(gene %in% gene_order) %>%
        pivot_longer(cols=c(n_up, n_dn), names_to="Direction", values_to="# case") %>% 
        mutate(`# case`=ifelse(Direction=="n_up", `# case`, -`# case`)) %>%
        mutate(Direction = ifelse(Direction=="n_up", "Increase in center", "Increase in periphery")) %>%
        mutate(gene = factor(gene, levels=rev(gene_order)))
    p = ggplot(all_cor_df_plot, aes(x=gene, y=`# case`, fill=cancer_type)) +
        geom_bar(stat="identity", position="stack", alpha=1) +
        theme_bw_nogrid +
#	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color="black",size=1)) +
	theme(axis.text = element_text(color="black",size=5))+
        coord_flip() +
	scale_fill_manual(values=cancer.color)+ 
	scale_y_continuous(labels=abs, breaks=scales::pretty_breaks()) +
        labs(title="Shared Layer DEGs") +
	geom_vline(xintercept=0)
    ggsave(str_interp("${output_dir}/${setting}_sample_in_case/${comp_type}/Layer_frac_case_rho_dir_count_${cancer}.pdf"), plot=p, w=2.5, h=3.5)
}
