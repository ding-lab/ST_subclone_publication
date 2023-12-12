# Run inferCNV on one ST sample
# 2023/06/02

library(infercnv)
library(Seurat)
library(reticulate)
library(googlesheets4)
library(tidyverse)
library(optparse)
# use_condaenv("inferCNV", required=TRUE)
reticulate::py_config()

sample_id = "HT206B1-S1Fc1U2Z1B1"
label_column = "Filtered_tumor_regions"
infercnv_mode = "sample"
output_prefix = ""
window_length = 151
sim_method = "meanvar"
num_normal = 200
pnorm_bayes = 0.3

### Parameters ------------------------------------------------------
option_list = list(
	  make_option(c("-s", "--sample_id"),
		type="character",
		default=NULL,
		help="Sample name for ST",
		metavar="character"),
	  make_option(c("-m", "--infercnv_mode"),
		type="character",
		default="subcluster",
		help="InferCNV analysis_mode sample|subcluster|cell",
		metavar="character"),
	  make_option(c("-c", "--label_column"),
		type="character",
		default="tumor_normal",
		help="InferCNV cell_type_file column name",
		metavar="character"),
	  make_option(c("-w", "--window_length"),
		type="integer",
		default="151",
		help="InferCNV smoothing window length",
		metavar="integer"),
	  make_option(c("-d", "--sim_method"),
		type="character",
		default="meanvar",
		help="InferCNV hspike simulation method",
		metavar="character"),
	  make_option(c("-n", "--num_normal"),
		type="integer",
		default="-1",
		help="Number normal spots to give input to inferCNV",
		metavar="integer"),
	  make_option(c("-p", "--pnorm_bayes"),
		type="numeric",
		default="0.5",
		help="Bayes posterior probability cutoff for inferCNV",
		metavar="numeric"),
	  make_option(c("-r", "--leiden_resolution"),
		type="numeric",
		default="0.05",
		help="Leiden resolution for clustering",
		metavar="numeric"),
	  make_option("--output_prefix",
	    type="character",
		default="",
		help="Output prefix",
		metavar="character")
	  )

	opt_parser = OptionParser(option_list=option_list);
	opt = parse_args(opt_parser)

sample_id=opt$sample_id
infercnv_mode=opt$infercnv_mode
label_column=opt$label_column
window_length=opt$window_length
sim_method=opt$sim_method
num_normal=opt$num_normal
pnorm_bayes=opt$pnorm_bayes
output_prefix=opt$output_prefix

short_label_column = case_when(label_column=="Filtered_tumor_regions"~"region",
			       label_column=="tumor_normal"~"TN",
			       TRUE~NA_character_)
normal_label = case_when(label_column=="Filtered_tumor_regions"~"Ref",
			 label_column=="tumor_normal"~"Normal",
			 TRUE~NA_character_)

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")

if (infercnv_mode == "sample") {
    output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_5-InferCNV_${output_prefix}w${window_length}_d${sim_method}_n${num_normal}")
} else if (infercnv_mode == "subcluster") {
    output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_5-InferCNV_${output_prefix}subcluster_w${window_length}_d${sim_method}_n${num_normal}")
} else if (infercnv_mode == "cell") {
    output_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_5-InferCNV_${output_prefix}cell_w${window_length}_d${sim_method}_n${num_normal}")
}

sample_outdir = str_interp("${output_dir}/${sample_id}/")
if (!dir.exists(sample_outdir)) {dir.create(sample_outdir)} 

# Don't rerun sample already done
final_mat_filename = str_interp("${sample_outdir}/infercnv.observations.txt")
stopifnot(!file.exists(final_mat_filename))

print(paste0("[inferCNV] Start running for sample: ", sample_id))

# Read metadata from annotation sheet
obj <- read_st_obj(sample_id)
st_meta = read_subclone_tbl(sample_id)
obj = AddMetaData(obj, st_meta[Cells(obj),])
obj$tumor_normal = ifelse(obj$Filtered_tumor_regions==0, "Normal", "Tumor")
if (label_column == "Filtered_tumor_regions") {
  subclone_color = create_subclone_color(obj$Filtered_tumor_regions)
  subclone_color["0"] = "grey75" 
  subclone_color["Ref"] = NA
} else if (label_column == "tumor_normal") {
  subclone_color = c("Tumor"= as.character(bound.color["Tumor"]), "Normal"=as.character(bound.color["TME"]))
}
Idents(obj) <- label_column
table(Idents(obj))

counts_matrix = GetAssayData(object = obj, slot = "counts")
gene_order = read.delim(file = "/diskmnt/Projects/Users/cliu/metnet_mCRC/scripts/InferCNV/gencode_v32_gene_name.txt", header=FALSE, stringsAsFactors = FALSE, sep="\t") 

cell_type_file = obj@meta.data %>% select(!!label_column)
gene_order_file = gene_order %>% column_to_rownames("V1")
Idents(obj) = label_column  

##### select normal spots with low estimated purity #####
if (num_normal > 0) {
    threshold <- sort(obj$Purity[obj$Filtered_tumor_regions==0])[num_normal]
    if (is.na(threshold)) {
	threshold = sort(obj$Purity)[num_normal]
        obj$Filtered_tumor_regions = ifelse(obj$Purity<=threshold, "Ref", obj$Filtered_tumor_regions)
    }
    obj$Filtered_tumor_regions = case_when(
        obj$Purity <= threshold & obj$Filtered_tumor_regions==0 ~ "Ref",
	TRUE ~ obj$Filtered_tumor_regions
    )
    cell_type_file_filtered = obj@meta.data %>% select(!!label_column)

	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
					    annotations_file=cell_type_file_filtered,
					    gene_order_file=gene_order_file,
					    ref_group_names=as.vector(normal_label))
} else {
	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
										annotations_file=cell_type_file,
										gene_order_file=gene_order_file,
										ref_group_names=as.vector(normal_label))
}
table(obj@meta.data[,label_column])

if (!dir.exists(str_interp("${output_dir}/infercnv_label"))) {dir.create(str_interp("${output_dir}/infercnv_label"))} 
p_HE = SpatialPlot(obj, group.by=label_column, alpha=0, crop=F) + 
	NoLegend() + labs(title=sample_id)
p = SpatialPlot(obj, group.by=label_column, stroke=0, label=T, label.size=2, label.box=F, crop=F, alpha=0.8) + 
	scale_fill_manual(values=subclone_color, na.value=NA) + NoLegend() +
 	labs(title="InferCNV label")
ggsave(str_interp("${output_dir}/infercnv_label/infercnv_label_${sample_id}.pdf"), plot=p_HE+p, w=6, h=3.5)

infercnv_obj = infercnv::run(infercnv_obj, 
			     window_length = window_length,
                             analysis_mode = infercnv_mode, 
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
			     out_dir=str_interp("${output_dir}/${sample_id}/"), # dir is auto-created for storing outputs 
                             cluster_by_groups=T,   # cluster
			     sim_method=sim_method,
                             plot_steps=T,
                             denoise=T,
                             HMM=T,
			     BayesMaxPNormal = pnorm_bayes, 
                             mask_nonDE_genes = T,
                             resume_mode=F,
                             num_threads = 8)

print(paste0("[inferCNV] Finished running for sample: ", sample_id))
