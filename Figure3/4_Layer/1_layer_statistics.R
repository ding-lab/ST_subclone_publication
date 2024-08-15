# Layer statistics for each Visium ST sample
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
library(patchwork)
library(ggpmisc)

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")
selected_tbl = read_tracking_sheet()

sample_id = "HT260C1-Th1K1Fc2U1Z1Bs1"
output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/40-Layer_DEG/40_1-Layer_DEG/")

### Set up sample to be processed -----------------------------------
case_id_list = c("HT260C1","HT397B1","HT268B1","HT270P1","HT112C1")
sample_id_list = selected_tbl %>% 
	filter(Study_cohort=="Discovery", case_id %in% case_id_list) %>%
	pull(Sample) %>% sort
sample_id_list = selected_tbl %>% filter(Study_cohort=="Discovery") %>% 
	pull(Sample) %>% sort

### Read clone df
clone_df = data.frame()
clone_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_13-InferCNV_settings/prob_jaccard_similarity_case/"
for (case_id in case_id_list) {
    clone_df = bind_rows(clone_df, 
			 read_tsv(str_interp("${clone_dir}/Clone_manualK_${case_id}.tsv"))
    )
}
clone_df$Filtered_tumor_regions = as.character(clone_df$Filtered_tumor_regions)

### Read all ST obj -------------------------------------------------
st_list = list()
p_layer_list = list()
p_big_layer_list = list()
p_big_fc3_list = list()
subclone_tbl_all = data.frame()
for (sample_id in sample_id_list) {
    st = read_st_obj(sample_id)
    subclone_tbl = read_subclone_tbl(sample_id) %>% 
	mutate(FC5 = ifelse(FC5=="unknown", NA_real_,as.numeric(FC5))) %>% 
	mutate_at(vars(starts_with("FC")), function(x){ifelse(x==0, NA_real_,x)})
    subclone_tbl_all = bind_rows(subclone_tbl_all, subclone_tbl %>% mutate(sample_id=sample_id))
    st = AddMetaData(st, subclone_tbl[Cells(st), ])
    p_HE = SpatialPlot(st, alpha=0, crop=F) + NoLegend() + labs(title=sample_id) 
    p = SpatialFeaturePlot(st, features=c("FC1","FC3","FC5"), crop=F, stroke=0, ncol=1) &
	    scale_fill_gradientn(colors=rev(sequential_hcl(palette="Inferno", 10)), na.value=NA)
    p_sample = wrap_plots(p_HE, p, ncol=1, heights=c(1,4))
    st_list[[sample_id]] = st
    p_layer_list[[sample_id]] = p_sample

    region_big = subclone_tbl %>% 
        filter(Filtered_tumor_regions!="0") %>%
        group_by(Filtered_tumor_regions) %>% 
        summarise(n_spot=n(),
	      FC1_max=max(FC1, na.rm=T), FC3_max=max(FC3, na.rm=T),
	      FC5_max=max(FC5, na.rm=T)) %>% 
        filter(n_spot>=50, (FC1_max>=3|FC3_max>=3|FC5_max>=3)) %>% 
        pull(Filtered_tumor_regions) %>% as.character

    st_tumor_big = subset(st, Filtered_tumor_regions %in% region_big)
#    st_tumor_big$FC3_max = st_tumor_big@meta.data %>% group_by(Filtered_tumor_regions) %>% 
#	    mutate(FC3_max = max(FC3, na.rm=T)) %>% pull(FC3_max)
#    st_tumor_big$FC3_frac = (st_tumor_big$FC3-1) / (st_tumor_big$FC3_max-1)

    p_big_fc3 = SpatialFeaturePlot(st_tumor_big, features=c("FC3","FC3_frac"), crop=F, stroke=0, ncol=1) &
	    scale_fill_gradientn(colors=rev(sequential_hcl(palette="Inferno", 10)), na.value=NA)
    p_big_fc3_list[[sample_id]] = wrap_plots(p_HE, p_big_fc3, ncol=1, heights=c(1,4/3*2))
    
    p_big = SpatialFeaturePlot(st_tumor_big, features=c("FC1", "FC3","FC5"), crop=F, stroke=0, ncol=1) &
	    scale_fill_gradientn(colors=rev(sequential_hcl(palette="Inferno", 10)), na.value=NA)
    p_big_layer_list[[sample_id]] = wrap_plots(p_HE, p_big, ncol=1, heights=c(1,4))
}

# pdf 
ggsave(str_interp("${output_dir}/Layers_big_fc3.png"), plot=wrap_plots(p_big_fc3_list, ncol=13), 
       w=3*length(sample_id_list)/3, h=30, dpi=72)

ggsave(str_interp("${output_dir}/Layers_big.png"), plot=wrap_plots(p_big_layer_list, ncol=13), 
       w=3*length(sample_id_list)/3, h=45, dpi=72)

ggsave(str_interp("${output_dir}/Layers.png"), plot=wrap_plots(p_layer_list, ncol=13), 
       w=3*length(sample_id_list)/3, h=45, dpi=72)

### FC3 max vs subclone size
sample_id_list = selected_tbl %>% pull(Sample) %>% sort
sample_id_list = selected_tbl %>% filter(Study_cohort=="Discovery") %>% 
	pull(Sample) %>% sort
subclone_tbl_all = data.frame()
for (sample_id in sample_id_list) {
    print(sample_id)
    subclone_tbl = read_subclone_tbl(sample_id) %>% 
	mutate(FC5 = ifelse(FC5=="unknown", NA_real_,as.numeric(FC5))) %>% 
	mutate_at(vars(starts_with("FC")), function(x){ifelse(x==0, NA_real_,x)})
    subclone_tbl_all = bind_rows(subclone_tbl_all, subclone_tbl %>% mutate(sample_id=sample_id))
}

df_toplot = subclone_tbl_all %>% 
    filter(Filtered_tumor_regions!="0") %>%
    left_join(selected_tbl %>% distinct(Sample, Study_cohort), by=c("sample_id"="Sample")) %>% 
    distinct(sample_id, Study_cohort, Filtered_tumor_regions) %>% 
    group_by(sample_id, Study_cohort) %>% 
    summarise(n_region=n())
p = ggplot(df_toplot, aes(x=sample_id, y=n_region)) +
        facet_wrap(~Study_cohort,scale="free_x") +
	geom_point() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))
ggsave(str_interp("${output_dir}/Region_summary.pdf"), plot=p, w=12, h=5)

layer_summary_df = subclone_tbl_all %>% 
    filter(Filtered_tumor_regions!="0") %>%
    group_by(sample_id, Filtered_tumor_regions) %>% 
    summarise(n_spot=n(), FC5_max=max(FC5, na.rm=T)) %>% 
#    left_join(clone_df %>% distinct(sample_id, clone_id, Filtered_tumor_regions)) %>% 
    left_join(selected_tbl %>% distinct(Sample, cancer_type, sample_type_tumor), by=c("sample_id"="Sample")) %>% 
    mutate(big_region = ifelse(n_spot>=50 & FC5_max>=3, "Big", "Small"))
p = ggplot(layer_summary_df, aes(x=n_spot, y=FC5_max, size=n_spot)) +
    geom_point(aes(alpha=big_region), shape=21) +
    facet_wrap(~sample_id, nrow=8, scales="free") +
    theme_bw_nogrid +
    scale_alpha_manual(values=c("Big"=1,"Small"=0.2)) +
    scale_y_continuous(breaks= scales::pretty_breaks()) +
    labs(x="Region size", size="Region size", y="Layer max")
ggsave(str_interp("${output_dir}/Layers_summary_bysample.pdf"), plot=p, w=10, h=10)

circle_ref = data.frame(y=c(1:12)) %>% 
#    mutate(x=pi*y**2) %>% 
    mutate(x=1+6*(y-1), ref="hexagon")
line_ref = data.frame(y=c(1:12)) %>% 
    mutate(x=y, ref="line") 

ref_df = bind_rows(circle_ref, line_ref) %>% 
    mutate(cancer_type=list(c("CRC","BRCA","PDAC"))) %>% 
    unnest(cancer_type)

max_layer = layer_summary_df %>% group_by(cancer_type) %>% summarise(y_max=max(FC5_max,na.rm=T))

p = ggplot(layer_summary_df, aes(x=n_spot, y=FC5_max)) +
    geom_jitter(aes(color=cancer_type, group=big_region, shape=sample_type_tumor),  
		width=0.5, height=0.5, size=1) +
    geom_line(data=circle_ref, aes(x=x,y=y), linetype="dashed") +
    facet_wrap(~cancer_type, nrow=1, scales="free") +
    theme_bw_nogrid +
    scale_color_manual(values=cancer.color, aesthetics=c("fill","color")) +
#    scale_size(limits=c(1,5)) +
#    scale_alpha_manual(values=c("Big"=1,"Small"=0.2)) +
    scale_shape_manual(values=c("Primary"=16, "Metastasis"=8)) +
    scale_y_continuous(breaks= scales::pretty_breaks()) +
    theme(text=element_text(size=16), legend.position="bottom") +
    labs(x="Region size (# of spots)", size="Region size", y="Layer max", fill="Cancer", color="Cancer", shape="Sample")
ggsave(str_interp("${output_dir}/Layers_summary.pdf"), plot=p, w=8, h=3.5)


