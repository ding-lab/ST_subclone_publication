# Identify CCI differentially enriched in the tumor boudnary regions from COMMOT outputs
# Jingxian Clara Liu
# 2023/10/08

library(Matrix)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggrepel)
suppressMessages(library(optparse))

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")

input_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/45-commotAnalysis/CellChat/"
sample_id = "HT112C1-U1_ST_Bn1"

option_list = list(
  make_option(c("-s", "--sample_id"),
	type="character",
	default=NULL,
	help="Sample name for ST",
	metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
sample_id=opt$sample_id

subclone_tbl = read_subclone_tbl(sample_id)

sample_output_dir = str_interp("${input_dir}/../CellChat_summary/${sample_id}/")
if (!dir.exists(sample_output_dir)) {dir.create(sample_output_dir)}

# Read COMMOT summary
test_dir="/diskmnt/Projects/Users/cliu/pancan_ST/TME/COMMOT/"
sender_tbl = read_delim(str_interp("${input_dir}/${sample_id}/sum-sender_${sample_id}.csv")) %>% 
    rename(barcode="...1")
receiver_tbl = read_delim(str_interp("${input_dir}/${sample_id}/sum-receiver_${sample_id}.csv")) %>% 
    rename(barcode="...1") 
barcode_anno_tbl = subclone_tbl %>% 
    mutate(spot_type = case_when(
	!is.na(Tumor_boundary) | !is.na(TME_boundary) ~ "Boundary",
	Filtered_tumor_regions != "0" ~ "Tumor",
	Filtered_tumor_regions == "0" ~ "TME"
    )) %>% 
    select(barcode, spot_type)
pathway_all_tbl = bind_rows(
    sender_tbl %>% pivot_longer(cols=!barcode, names_to="pathway", values_to="x") %>% mutate(direction="sender"),
    receiver_tbl %>% pivot_longer(cols=!barcode, names_to="pathway", values_to="x") %>% mutate(direction="receiver")
    ) %>%
    left_join(barcode_anno_tbl, by="barcode")  
pathway_summary_tbl = pathway_all_tbl %>%  
#    group_by(pathway, spot_type) %>%
    group_by(pathway) %>%
    summarise(
	max_signal = max(x),
	mean_signal = mean(x),
	median_signal = median(x),
	quant90_signal = quantile(x, 0.9, na.rm=T)
    ) %>% 
    ungroup %>% 
    filter(!grepl("total-total", pathway))
pathway_summary_tbl %>% slice_max(n=50, order_by=quant90_signal) %>% as.data.frame

####### Boundary enriched COMMOT signal #########
if (!file.exists(str_interp("${sample_output_dir}/Volcano_COMMOT_${sample_id}.pdf"))) {
pathway_de_p = pathway_all_tbl %>%
    filter(!grepl("total-total", pathway)) %>%
    mutate(boundary = factor(ifelse(spot_type=="Boundary", "B","NB"), levels=c("B","NB"))) %>% 
    group_by(pathway, direction) %>%
    summarise(pval = wilcox.test(x~boundary, paired=FALSE)$p.value) %>% 
    mutate(p.val.adj = p.adjust(pval, method="fdr")) %>% 
    ungroup 
# pathway_median = pathway_summary_tbl %>% 
#    pivot_wider(id_cols=pathway, names_from="spot_type", values_from="median_signal")
pathway_median = pathway_all_tbl %>% 
    filter(!grepl("total-total", pathway)) %>%
    mutate(boundary = factor(ifelse(spot_type=="Boundary", "B","NB"), levels=c("B","NB"))) %>%
    group_by(pathway, direction, boundary) %>%
    summarise(median_signal = median(x)) %>%
    pivot_wider(id_cols=c(pathway,direction), names_from="boundary", values_from="median_signal")
pathway_de_tbl = pathway_de_p %>% 
    left_join(pathway_median, by=c("pathway","direction")) %>%
    arrange(p.val.adj, pval) 
pathway_de_tbl %>% head
pathway_de_tbl %>% arrange(desc(B-NB))

pathway_de_tbl %>% 
    write_tsv(str_interp("${sample_output_dir}/DE_COMMOT_B_NB_${sample_id}.tsv"))

df_toplot = pathway_de_tbl %>%
#    mutate(rank_stats = abs(B-NB)*(-log10(p.val.adj))) %>%
    mutate(rank_stats = abs(B-NB)) %>% 
    arrange(desc(rank_stats)) %>%
    group_by(B>NB) %>% 
    mutate(label = ifelse(row_number()<=15 & p.val.adj<0.05 & abs(B-NB)>0.2, pathway, NA_character_)) 
p = ggplot(df_toplot, aes(x=B-NB, y=-log10(p.val.adj))) +
	geom_point(aes(color=B-NB)) +
	geom_hline(yintercept = -log10(0.05), linetype="dashed") +
	geom_vline(xintercept = 0) +
	theme_classic() +
	theme(axis.text=element_text(color="black",size=16),
	      axis.title=element_text(color="black",size=16)) +
	ggrepel::geom_text_repel(aes(label=label), size=2, force_pull=10, 
                    segment.size = 0.1, min.segment.length = 0,
                    max.overlaps = Inf, ylim=c(-log10(0.05), NA)) +
	scale_color_gradient2(high="firebrick", low="dodgerblue", mid="grey75") +
	labs(title=sample_id)
ggsave(str_interp("${sample_output_dir}/Volcano_COMMOT_${sample_id}.pdf"), w=5, h=5, plot=p)
}

####### All COMMOT signal distribution #######
if (!file.exists(str_interp("${sample_output_dir}/Boxplot_QC_COMMOT_${sample_id}.pdf"))) {
df_toplot = pathway_summary_tbl %>% 
    pivot_longer(cols=ends_with("_signal"), names_to="stats", values_to="value") %>% 
    group_by(stats) %>% 
    arrange(desc(value)) %>%
    mutate(label=ifelse(row_number()<=10, pathway, NA_character_))
# plot_max = quantile(df_toplot$value, 0.99) 
p = ggplot(df_toplot, aes(x=stats, y=log2(value)))+
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(size=0.2, width=0.2, stroke=0) +
    geom_text(aes(label=label), size=1) + 
#    ylim(c(0, plot_max)) +
    theme_classic() +
    theme(axis.text=element_text(color="black",size=12),
	  axis.title=element_text(color="black",size=12)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(str_interp("${sample_output_dir}/Boxplot_QC_COMMOT_${sample_id}.pdf"), w=4, h=5, plot=p)
}

############# Time consuming!! #############
if (F) {
# Read COMMOT output
sample_dir = str_interp("${input_dir}/${sample_id}/")
mat_list = list()
lr_tbl_long = data.frame()
for (mtx_path in list.files(sample_dir)) {
    pair = gsub(".*_commot-cellchat-(.*)\\.mtx","\\1",mtx_path)
    print(pair)
    abs_mtx_path = str_interp("${sample_dir}/${mtx_path}")
    mat = readMM(abs_mtx_path)
    mat_list[[pair]] = mat
    lr_tbl_long_sample = mat %>% 
	summary %>% as.data.frame %>%
	mutate(pair = pair)
    lr_tbl_long = bind_rows(lr_tbl_long, lr_tbl_long_sample)
}

# Annotate LR summary table with boundary type
boundary_i = subclone_tbl %>% 
    mutate(i=row_number()) %>%
    filter(!is.na(Tumor_boundary) | !is.na(TME_boundary)) %>%
    pull(i)
tumor_i = which(subclone_tbl$Filtered_tumor_regions!="0")
tme_i = which(subclone_tbl$Filtered_tumor_regions=="0")
lr_tbl_anno = lr_tbl_long %>%
    mutate(spot_spot_type = case_when(
	(i %in% boundary_i) & (j %in% boundary_i) ~ "B-B",
	(i %in% tumor_i) & (j %in% tumor_i) ~ "NB-NB_tumor",
	(i %in% tme_i) & (j %in% tme_i) ~ "NB-NB_tme",
	TRUE ~ "B-NB"
    ))  
total_tbl_anno = lr_tbl_anno %>% filter(pair=="total-total")
lr_tbl_anno = lr_tbl_anno %>% filter(pair!="total-total")
lr_tbl_summary = lr_tbl_anno %>%
#    group_by(pair, spot_spot_type) %>%
    group_by(pair) %>% 
    summarise(
	max_signal = max(x),
	mean_signal = mean(x),
	median_signal = median(x)
    ) %>% 
    ungroup
lr_tbl_summary %>% slice_max(n=20, order_by=max_signal)

# Signal level total-total
p1 = ggplot(lr_tbl_anno, aes(x=x, fill=spot_spot_type)) +
    geom_histogram(position="dodge") +
    scale_y_continuous(trans="log1p") +
    theme_bw() +
    labs(title="All pairs")
p2 = ggplot(total_tbl_anno, aes(x=x, fill=spot_spot_type)) +
    geom_histogram(position="dodge") +
    scale_y_continuous(trans="log1p") +
    theme_bw() +
    labs(title="total-total")
ggsave(str_interp("${sample_output_dir}/Histogram_SpotBoundary_${sample_id}.pdf"), 
       plot=p1/p2, w=6, h=6)
}
############################################
