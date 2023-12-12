# RCTD subclone TME composition summary at the tumor boundary regions
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

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")
selected_tbl = read_tracking_sheet()
input_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/13-RCTD/13_2-subclone_summary/best_version/"
output_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/13-RCTD/13_3-cohort/Sample_ref/ByLayer/"

### Set up sample to be processed -----------------------------------
sample_id_list = list.dirs(input_dir, full.names=F) %>% setdiff(c(""))
sample_id_list = list.files(input_dir, full.names=F, pattern="*_RCTD_weights.tsv") %>% gsub("_RCTD_weights.tsv","",.)

immune_celltypes = c(cell_types_myeloid, cell_types_lymphoid)
clean_id_mapping = selected_tbl$SampleID_Clean_U %>% `names<-`(selected_tbl$Sample)

region_size_df = read_tsv(REGION_SUMMARY_PTH) 

clone_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/"
clone_df = read_tsv(str_interp("${clone_dir}/Consensus_genetic_clones_2023-10-05.tsv")) %>% 
	mutate(clone_id = LETTERS[as.numeric(gsub("clone_","",genetic_clone))])
clone_df$sample_id = clean_id_mapping[clone_df$sample_id]

### Read all RCTD composition --------------------------------------------
df_bound_raw = data.frame()
sample_id = "HT260C1-Th1K1Fc2U1Z1Bs1"
sample_id = "HT308B1-S1H4Fc2U1Z1Bs1"
for (sample_id in sample_id_list) {
  file_name = str_interp("${input_dir}/${sample_id}_RCTD_weights.tsv")
  if (!file.exists(file_name)) {next}
  subclone_tbl = read_tsv(file_name)

  for (region in sort(unique(subclone_tbl$t)))

  ## Read and process RCTD summary for inner tumor
  df_tumor_sample = read_tsv(file_name) %>%
	  mutate(region_type=case_when(
			!is.na(boundary_type) ~ NA_character_,
			Filtered_tumor_regions>0 ~ "Tumor",
			Filtered_tumor_regions==0 ~ "TME",
			TRUE ~ NA_character_)) %>%
  	  select(-boundary_type) %>%
	  filter(!is.na(region_type)) %>%
	  column_to_rownames("barcode") %>%
	  pivot_longer(cols=!c("region_type", "Filtered_tumor_regions"), names_to="celltype", values_to="RCTD_weight") %>%
	  mutate(celltype = gsub("RCTD_","",celltype)) %>%
	  group_by(Filtered_tumor_regions, region_type, celltype) %>%
	  summarise(RCTD_weight_sum = sum(RCTD_weight, na.rm=T)) %>%
	  group_by(Filtered_tumor_regions, region_type) %>% 
	  mutate(RCTD_weight_sum = RCTD_weight_sum/sum(RCTD_weight_sum, na.rm=T)) %>%
	  mutate(sample_id = sample_id) 

  ## Read RCTD bound summary
  file_name = str_interp("${output_dir}/${sample_id}/${sample_id}_RCTD_weights_bound_summary.tsv")
  if (!grepl("HT308B1", sample_id)) {file_name=gsub("Sample_ref_test","Sample_ref",file_name)}
  if (!file.exists(file_name)) {next}
  
  df_bound_sample = read_tsv(file_name) %>%
	  mutate(sample_id = sample_id) 
  df_bound_raw = bind_rows(df_bound_raw, df_tumor_sample) %>% bind_rows(df_bound_sample)
}

# Clean up RCTD bound summary
df_bound_raw = df_bound_raw %>%
  mutate(celltype=gsub("s$","",celltype)) %>%
  mutate(celltype = case_when(
	celltype %in% names(cell_type_mapping) ~ cell_type_mapping[celltype],
	TRUE ~ celltype
  )) %>%  
  mutate(sample_id = clean_id_mapping[sample_id]) %>%
  rename(RCTD_weight = "RCTD_weight_sum")

df_bound = df_bound_raw 
df_fill0 = df_bound %>% expand(nesting(Filtered_tumor_regions, sample_id, region_type, celltype)) 
df_bound = df_bound %>% right_join(df_fill0) %>%
	mutate(RCTD_weight=replace_na(RCTD_weight, 0)) %>%
	mutate(sample_tumor_regions = paste0(sample_id, "_", Filtered_tumor_regions)) 

sample_id_order_df = selected_tbl %>% filter(SampleID_Clean_U%in%df_bound$sample_id) %>% 
	arrange(cancer_type, SampleID_Clean_U) 
sample_id_order = sample_id_order_df$SampleID_Clean_U
df_bound$sample_id = factor(df_bound$sample_id, levels=sample_id_order)

cell_type_names = unique(df_bound$celltype)

sample.colors = cancer.color[sample_id_order_df$cancer_type] %>% `names<-`(NULL)

### RCTD scatterplot boundary type
df_toplot = df_bound %>%
	mutate(region_type = factor(region_type, levels=c("Tumor",boundary_level,"TME"))) %>%
	group_by(region_type, celltype, sample_id) %>%
	summarise(RCTD_weight_sum = sum(RCTD_weight, na.rm=T)) %>%
	ungroup %>%
	group_by(region_type, sample_id) %>%
	mutate(RCTD_weight_sum = RCTD_weight_sum/sum(RCTD_weight_sum, na.rm=T)) %>% 
	mutate(sample_celltype = paste0(sample_id, "_", celltype))

p <- ggplot(df_toplot, aes(x=region_type, y=RCTD_weight_sum, color=celltype, group=sample_celltype)) +
	geom_line() +
	geom_point(aes(size=RCTD_weight_sum)) +
	facet_grid(~ sample_id, scale="free", space="free") +
	scale_size(range=c(0.5,2)) +
        scale_color_manual(values=celltype.color, breaks=names(celltype.color)) +
        scale_y_continuous(labels=abs) +
        theme_bw_nogrid +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6)) +
	guides(size=guide_legend(nrow=1)) +
        labs(title=str_interp("Cell type composition: trend across boundary"), x="Boundary type", y="Fraction",
	     color="Cell type", size="RCTD weight")
ggsave(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary.pdf"), 
       plot=color_strip(p, sample.colors), w=18, h=6)

### RCTD scatterplot boundary type - split by celltype
df_toplot = df_bound %>%
	left_join(clone_df) %>%
	mutate(sample_clone = paste0(sample_id,"_",clone_id)) %>% 
	mutate(region_type = factor(region_type, levels=c("Tumor",boundary_level,"TME"))) %>%
	group_by(region_type, celltype, sample_id, sample_clone) %>%
	summarise(RCTD_weight_sum = sum(RCTD_weight, na.rm=T)) %>%
	ungroup %>%
	group_by(region_type, sample_id, sample_clone) %>%
	mutate(RCTD_weight_sum = RCTD_weight_sum/sum(RCTD_weight_sum, na.rm=T)) %>% 
	mutate(sample_celltype = paste0(sample_id, "_", celltype))
df_toplot = df_bound %>%
	left_join(clone_df) %>%
	mutate(sample_clone = paste0(sample_id,"_",clone_id)) %>% 
	mutate(region_type = factor(region_type, levels=c("Tumor",boundary_level,"TME"))) %>%
	mutate(RCTD_weight_sum=RCTD_weight)
shared_celltypes = df_toplot %>% ungroup %>% 
	group_by(celltype) %>% mutate(max_frac=max(RCTD_weight_sum)) %>% ungroup %>% 
	distinct(sample_id, celltype, max_frac) %>%
	group_by(celltype, max_frac) %>%
	summarise(n_sample = n()) %>% 
	filter(n_sample>5, max_frac>0.01) %>% pull(celltype)
df_toplot = df_toplot %>% filter(celltype %in% shared_celltypes)

p <- ggplot(df_toplot, aes(x=region_type, y=RCTD_weight_sum, color=celltype)) +
	geom_boxplot(outlier.shape=NA) +
#	geom_line(aes(group=sample_id)) +
#	geom_point(aes(size=RCTD_weight_sum)) +
	facet_wrap(~ celltype, scale="free", nrow=4) +
	scale_size(range=c(0.5,2)) +
        scale_color_manual(values=celltype.color, breaks=names(celltype.color)) +
        scale_y_continuous(labels=abs) +
        theme_bw_nogrid +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6), legend.position="bottom") +
#	guides(size=guide_legend(nrow=1)) +
        labs(title=str_interp("Cell type composition: trend across boundary"), x="Boundary type", y="Fraction",
	     color="Cell type", size="RCTD weight")
ggsave(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_split_bycelltype.pdf"), 
       plot=p, w=8, h=12)

celltype="Macrophage"
celltype="T"
p <- ggplot(df_toplot %>% filter(celltype==!!celltype), 
	    aes(x=region_type, y=RCTD_weight_sum, color=celltype, group=sample_clone)) +
	geom_boxplot(outlier.shape=NA) +
#	geom_line() +
#	geom_point(aes(size=RCTD_weight_sum)) +
	facet_wrap(~ sample_id, nrow=4, scale="free") +
	scale_size(range=c(0.5,2)) +
        scale_color_manual(values=celltype.color, breaks=names(celltype.color)) +
        scale_y_continuous(labels=abs) +
        theme_bw_nogrid +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6), legend.position="bottom") +
#	guides(size=guide_legend(nrow=1)) +
        labs(title=str_interp("Cell type composition: trend across boundary"), x="Boundary type", y="Fraction",
	     color="Cell type", size="RCTD weight")
ggsave(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_split_${celltype}.pdf"), 
       plot=p, w=12, h=8)

### RCTD scatterplot boundary type -- per subclone
df_toplot = df_bound %>%
	left_join(clone_df) %>%
	mutate(celltype=ifelse(celltype%in%names(shared_celltype), shared_celltype[celltype], celltype)) %>%
	filter(celltype %in% c(shared_celltype,"Tumor")) %>%
	filter(!region_type%in%c("TME")) %>%
	mutate(region_type = factor(region_type, levels=c("Tumor",boundary_level,"TME"))) %>%
	group_by(region_type, celltype, sample_id, clone_id) %>%
	summarise(RCTD_weight_sum = sum(RCTD_weight, na.rm=T)) %>%
	ungroup %>%
	group_by(region_type, sample_id, clone_id) %>%
	mutate(RCTD_weight_sum = RCTD_weight_sum/sum(RCTD_weight_sum, na.rm=T)) %>% 
	mutate(sample_celltype = paste0(sample_id, "_", celltype), 
	       sample_region=paste0(sample_id,"_", clone_id)) 
celltype.color.toplot = celltype.color[c(unique(shared_celltype),"Tumor")] 

p <- ggplot(df_toplot, aes(x=region_type, y=RCTD_weight_sum, color=celltype, group=sample_region)) +
	geom_line(alpha=0.5) +
	geom_point(size=0.5, alpha=0.5) +
	facet_grid(celltype ~ sample_id, scale="free", space="free_x") +
	scale_size(range=c(0.5,2)) +
        scale_color_manual(values=celltype.color.toplot, breaks=names(celltype.color.toplot)) +
        scale_y_continuous(labels=abs) +
        theme_bw_nogrid +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6)) +
	guides(size=guide_legend(nrow=1)) +
        labs(title=str_interp("Cell type composition: trend across boundary"), x="Boundary type", y="Fraction",
	     color="Cell type", size="RCTD weight")
ggsave(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_bysubclone.pdf"), 
       plot=color_strip(p, sample.colors), w=18, h=12)

trend_df = data.frame()
for (sample_id in unique(df_toplot$sample_id)) {
    print(sample_id) 
    sample_df = df_toplot %>% filter(sample_id==!!sample_id)
    for (celltype in unique(sample_df$celltype)) {
	for (region in unique(sample_df$Filtered_tumor_regions)) {
	    unit_df = sample_df %>% filter(celltype==!!celltype, Filtered_tumor_regions==region)
    	    rctd = unit_df %>% ungroup %>% 
		full_join(unit_df %>%ungroup %>% expand(region_type)) %>% 
		mutate(RCTD_weight_sum=replace_na(RCTD_weight_sum, 0)) %>% 
		arrange(region_type) %>% pull(RCTD_weight_sum)
	    rctd_trend = c("0"="0","-1"="-","1"="+")[as.character(sign(rctd[2:6]-rctd[1:5]))] %>% paste(collapse="")
	    trend_df = bind_rows(
		trend_df,
		data.frame(sample_id=sample_id, celltype=celltype, region=region, trend=rctd_trend)
	    )
	}
    }
}
# trend_df %>% write_tsv(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_bysubclone_bycelltype_trend.tsv"))

region_size_df = data.frame()
for (sample_id in unique(df_toplot$sample_id)) {
    raw_sample_id = selected_tbl %>% filter(SampleID_Clean_U==sample_id) %>% pull(Sample)
    subclone_tbl = read_subclone_tbl(raw_sample_id)
    region_size_df = bind_rows(region_size_df,
			       subclone_tbl %>% group_by(Filtered_tumor_regions) %>% summarise(n_spot=n()) %>% mutate(sample_id=sample_id))
}

df_toplot_big = df_toplot
df_toplot_big = df_toplot %>% ungroup %>%
	mutate(Filtered_tumor_regions=as.character(Filtered_tumor_regions)) %>%
	left_join(region_size_df) %>% 
	filter(n_spot>=10)

p <- ggplot(df_toplot_big, 
	aes(x=region_type, y=RCTD_weight_sum, color=celltype, group=sample_region)) +
	geom_line(alpha=0.5) +
#	geom_point(size=0.5, alpha=0.5) +
	facet_wrap(~celltype, scale="free", nrow=2) +
	scale_size(range=c(0.5,2)) +
        scale_color_manual(values=celltype.color.toplot, breaks=names(celltype.color.toplot)) +
        scale_y_continuous(labels=abs) +
        theme_bw_nogrid +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6)) +
	guides(size=guide_legend(nrow=1)) +
        labs(title=str_interp("Cell type composition: trend across boundary"), x="Boundary type", y="Fraction",
	     color="Cell type", size="RCTD weight")
ggsave(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_bysubclone_bycelltype.pdf"), 
       plot=p, w=8, h=6)

trend_df = read_tsv(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_bysubclone_bycelltype_trend.tsv"))
trend_df = trend_df %>% 
	mutate(trend = str_sub(trend,0,-2))
df_toplot_trend_big = df_toplot_big %>%
	left_join(trend_df %>% mutate(region=as.character(region)), by=c("celltype"="celltype","Filtered_tumor_regions"="region","sample_id"="sample_id"))
df_toplot_trend_big$trend %>% table %>% sort

p <-  ggplot(df_toplot_trend_big, aes(x=celltype, fill=celltype)) +
	geom_bar() +
        scale_fill_manual(values=celltype.color.toplot, breaks=names(celltype.color.toplot)) +
	theme_bw() +
	theme(axis.text.x = element_blank()) +
	facet_wrap(~trend, nrow=3)
ggsave(str_interp("${output_dir}/Barplot_tumor_subclones_RCTD_weights_Boundary_bysubclone_bycelltype_trend.pdf"), 
       plot=p, w=12, h=8)

trend_list = c("++++","+---","----", "+++-")
df_toplot_trend_big = df_toplot_trend_big %>% 
	mutate(trend_new = case_when(
	    trend %in% c("----") ~ "----",
	    trend %in% c("+---") ~ "+---",
	    trend %in% c("++++","0+++","00++","000+") ~ "[0+]*+",
	    trend %in% c("+---", "+--0","+-00") ~ "+-[-0]*",
	    trend %in% c("+++-") ~ "+++-",
	    TRUE ~ "Other"
	))

p <- ggplot(df_toplot_trend_big,
	aes(x=region_type, y=RCTD_weight_sum, color=celltype, group=sample_region)) +
	geom_line(alpha=0.2) +
#	geom_point(size=0.5, alpha=0.5) +
	facet_grid(trend_new~celltype, scales="free") +
	scale_size(range=c(0.5,2)) +
        scale_color_manual(values=celltype.color.toplot, breaks=names(celltype.color.toplot)) +
        scale_y_continuous(labels=abs) +
        theme_bw_nogrid +
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6)) +
	guides(size=guide_legend(nrow=1)) +
        labs(title=str_interp("Cell type composition: trend across boundary"), x="Boundary type", y="Fraction",
	     color="Cell type", size="RCTD weight")
ggsave(str_interp("${output_dir}/Tumor_subclones_RCTD_weights_Boundary_bysubclone_bycelltype_trend.pdf"), 
       plot=p, w=12, h=8)

