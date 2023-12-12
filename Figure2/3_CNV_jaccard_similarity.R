# Calculate modified Jaccard similarity between two genome-wide CNV profiles 
# 2023/01/26
# Jingxian Clara Liu

library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
library(colorspace)
library(circlize)
library(optparse)
library(infercnv)
library(doParallel)
library(foreach)
library(EnrichedHeatmap)
library(patchwork)

set.seed(1234)

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")

# Test parameters
setting = "w151_dmeanvar_n200"
case_id = "HT112C1"
case_id = "HT397B1"
case_id = "HT262B1"
MAX_K=10
linkage_method = "ward.D2"
output_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_17-FFPE_workflow_v2/"

hatchet_dir = "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/Cong/hatchet_wes_v6/"
BLOCK_INFERCNV_SAMPLES = read_lines(str_interp("${ROOT_DIR}/23-genomic/23_17-FFPE_workflow_v2/BLOCK_INFERCNV_SAMPLES.txt"))
WHITELIST_WES_SAMPLES = read_lines(str_interp("${ROOT_DIR}/23-genomic/23_17-FFPE_workflow_v2/WHITELIST_WES_SAMPLES.txt"))

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
sample_id_list = selected_tbl %>% filter(case_id==!!case_id) %>% pull(Sample) %>% sort  %>% 
	setdiff(BLOCK_INFERCNV_SAMPLES)
WES_sample_id = read_column_val(sample_id_list[1], "WES_sample_id")
if (!WES_sample_id %in% WHITELIST_WES_SAMPLES) {WES_sample_id="No_WES"}
calico_sample_id = read_column_val(sample_id_list[1], "calico_sample_id")
sn_sample_id = read_column_val(sample_id_list[1], "Paired_snRNA_id")
hatchet_run = grep(case_id, dir(hatchet_dir), value=T)

output_dir = str_interp("${output_dir}/${case_id}")
if (!dir.exists(output_dir)) {dir.create(output_dir)}

###################################################
######## CNV Jaccard Similarity Helper ############
###################################################
# Similarity metrics between two segment-level CNV 
## Both inputs as GRange object with numeric metadata column CNV 
CNV_profile_jaccard = function(A_gr, B_gr, mode="w"){
  A_gr_pos = A_gr[A_gr$CNV>1]                                          
  A_gr_neg = A_gr[A_gr$CNV<1]                                          
  B_gr_pos = B_gr[B_gr$CNV>1]                                          
  B_gr_neg = B_gr[B_gr$CNV<1]
  intersect_gr_list = list()
  
  if (length(A_gr_pos)>0 & length(B_gr_pos)>0) {intersect_gr_list = c(intersect_gr_list, intersect(A_gr_pos, B_gr_pos))}
  if (length(A_gr_neg)>0 & length(B_gr_neg)>0) {intersect_gr_list = c(intersect_gr_list, intersect(A_gr_neg, B_gr_neg))}
  intersect_gr = do.call("c", intersect_gr_list)

  if (length(A_gr[A_gr$CNV!=1])>0) {
	union_gr = reduce(c(A_gr[A_gr$CNV!=1], B_gr[B_gr$CNV<1]))
  } else {
	union_gr = NULL
  }

  if (is.null(intersect_gr) | length(union_gr)==0 | length(intersect_gr)==0) {
	ratio = 0
  } else {
	ratio = sum(width(intersect_gr)) / sum(width(union_gr))
  }

  return( ratio )
}

# helper for inner product of two probability vectors as string separated by "," 
# e.g. A="0.1,0.2,0.7", B="0.3,0.4,0.3", return 0.1*0.3+0.2*0.4+0.7*0.3
helper_multi = function(A,B) {
  A = str_split(A,",") %>% unlist %>% as.numeric
  B = str_split(B,",") %>% unlist %>% as.numeric
  return(sum(A*B))
}

# Probablisitic similarity metrics between two segment-level CNV 
## Both inputs as GRange object with numeric metadata column CNV and numeric probability metadata column prob
CNV_profile_jaccard_prob = function(A_gr, B_gr, mode="all") {
  if (mode=="pred") {default_prob=0} else if (mode=="all") {default_prob="0,0,0"}
  gr_combined_list = list()
  for (state in c("amp","del")) {
    A_gr_pos = A_gr[A_gr$CNV_sign==state]
    B_gr_pos = B_gr[B_gr$CNV_sign==state]

    if(length(A_gr_pos)==0)  {gr_pos=B_gr_pos} else {gr_pos = disjoin(c(A_gr_pos, B_gr_pos))}
    if(length(gr_pos)==0) {gr_combined_list[[state]] = gr_pos; next}
    gr_pos$CNV_sign_A = 'neu'
    gr_pos$CNV_sign_B = 'neu'
    gr_pos$prob_A = default_prob
    gr_pos$prob_B = default_prob
    gr_pos$CNV_sign_A[to(findOverlaps(A_gr_pos, gr_pos))] = A_gr_pos[from(findOverlaps(A_gr_pos, gr_pos))]$CNV_sign 
    gr_pos$CNV_sign_B[to(findOverlaps(B_gr_pos, gr_pos))] = B_gr_pos[from(findOverlaps(B_gr_pos, gr_pos))]$CNV_sign 
    gr_pos$prob_A[to(findOverlaps(A_gr_pos, gr_pos))] = A_gr_pos[from(findOverlaps(A_gr_pos, gr_pos))]$prob 
    gr_pos$prob_B[to(findOverlaps(B_gr_pos, gr_pos))] = B_gr_pos[from(findOverlaps(B_gr_pos, gr_pos))]$prob

    gr_combined_list[[state]] = gr_pos
  }
  if(length(gr_combined_list[["amp"]])==0) {
    gr_combined = c(gr_combined_list[["del"]], gr_combined_list[["amp"]]) 
  } else {
    gr_combined = c(gr_combined_list[["amp"]], gr_combined_list[["del"]])
  }
  
  if (length(gr_combined) == 0) {
    ratio = 0
  } else {
    w_union = sum(width(gr_combined))
    if (mode=="pred") {
      w_intersect = sum(width(gr_combined)*
			as.numeric(gr_combined$CNV_sign_A==gr_combined$CNV_sign_B)*
		        gr_combined$prob_A*gr_combined$prob_B)
    } else if (mode=="all") {
      intersect_df = as.data.frame(gr_combined) %>%
	rowwise %>% mutate(prob_prod = helper_multi(prob_A,prob_B)) 
      w_intersect = sum(intersect_df$width*
			as.numeric(intersect_df$CNV_sign_A==intersect_df$CNV_sign_B)*
			intersect_df$prob_prod)
    }
    ratio = w_intersect/w_union
  }
}

###################################################
################# Read InferCNV ###################
###################################################
gr_list = list()
gr_pred_list = list()
gr_all_list = list()

for (sample_id in sample_id_list) {
subclone_tbl = read_subclone_tbl(sample_id) 

# Read inferCNV per region seg-level output
infercnv_dir = str_interp("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_5-InferCNV_${setting}/${sample_id}/")
infercnv_path = str_interp("${infercnv_dir}/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.3.pred_cnv_regions.dat")
mcmc_path = str_interp("${infercnv_dir}/BayesNetOutput.HMMi6.hmm_mode-samples/MCMC_inferCNV_obj.rds")

mcmc = readRDS(mcmc_path)
cnv_means = sapply(mcmc@cnv_probabilities,function(i) colMeans(i))
colnames(cnv_means) <- sapply(mcmc@cell_gene, function(i) { as.character(i$cnv_regions) })
rownames(cnv_means) <- c(sprintf("State:%s",seq_len(6)))

infercnv_hmm = read_tsv(infercnv_path) %>%
	mutate(CNV = infercnv_cnv_mapping_fn(state)) %>%
	mutate(CNV_sign = case_when(CNV>1~"amp", CNV<1~"del",CNV==1~"neu")) %>%
	mutate(Filtered_tumor_regions = gsub("\\..*","", cell_group_name))  
# Make sure each cnv_name is unique
stopifnot(nrow(infercnv_hmm) == nrow(distinct(infercnv_hmm, cell_group_name, cnv_name)))

prob_df = cnv_means %>% as.data.frame %>%  
	rownames_to_column("state") %>% 
	pivot_longer(cols=!c(state), names_to="cnv_name", values_to="prob") %>% 
	mutate(CNV = infercnv_cnv_mapping_fn(gsub("State:","",state))) %>%
	mutate(CNV_sign = case_when(CNV>1~"amp", CNV<1~"del",CNV==1~"neu")) %>%
	group_by(CNV_sign, cnv_name) %>% 
	summarise(prob=sum(prob, na.rm=T)) %>% ungroup
prob_df_vec = prob_df %>% pivot_wider(names_from="CNV_sign", names_prefix="prob_", values_from="prob") %>%
	mutate(prob = paste(prob_amp, prob_del, prob_neu, sep=",")) %>% 
	select(cnv_name, prob)

for (region in tumor_region_this_sample) {
    hmm = infercnv_hmm %>% filter(Filtered_tumor_regions==region)
    hmm_gr = GRanges(seqnames = hmm$chr, ranges = IRanges(hmm$start, hmm$end),
		     CNV_sign = hmm$CNV_sign, CNV=hmm$CNV)
    gr_list[[paste0(sample_id, "__", region)]] = hmm_gr

    ## only prob of predicted state 
    prob_r_pred = hmm %>% select(cnv_name, CNV_sign) %>% left_join(prob_df, by=c("cnv_name","CNV_sign"))
    hmm_gr_pred = hmm_gr 
    mcols(hmm_gr_pred)$prob = prob_r_pred$prob
    gr_pred_list[[paste0(sample_id, "__", region)]] = hmm_gr_pred

    ## prob sum of all possible states
    prob_r_all = hmm %>% select(cnv_name) %>% left_join(prob_df_vec, by=c("cnv_name")) %>% 
	    column_to_rownames("cnv_name")
    hmm_gr_all = hmm_gr
    mcols(hmm_gr_all) = cbind(mcols(hmm_gr_all), DataFrame(prob_r_all))
    gr_all_list[[paste0(sample_id, "__", region)]] = hmm_gr_all
}
}

########### IF WES, Filter WES regions ############
if (WES_sample_id!="No_WES") { 
    input_dir <- "/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/23-genomic/23_4-WES/20230629_ST/"
    gene_pos_filename <- str_interp("${input_dir}/gencode.v34.annotation.gene_filtered.protein_coding.ensembl_ID_no_version.protein-coding_hgnc_filtered.duplicates_removed.ensembl_ID_removed.txt")
    wes_df = read_tsv(str_interp("${input_dir}/${WES_sample_id}.T.called.seg"), comment="@") %>% 
    	filter(CALL!="0") %>% 
    	mutate(MEAN_LOG2_COPY_RATIO = ifelse(CALL!="0",2**MEAN_LOG2_COPY_RATIO,1)) %>%
    	dplyr::rename(Chromosome="CONTIG", Start="START", End="END",
    			Segment_Mean="MEAN_LOG2_COPY_RATIO")
    wes_gr <-  GRanges(seqnames = wes_df$Chromosome, ranges = IRanges(wes_df$Start, wes_df$End), CNV=wes_df$CALL)

    tumor_region_this_case = names(gr_list)
    for (region in tumor_region_this_case) {
        gr_list[[region]] = gr_list[[region]][unique(from(findOverlaps(gr_list[[region]], wes_gr)))] 
        gr_pred_list[[region]] = gr_pred_list[[region]][unique(from(findOverlaps(gr_pred_list[[region]], wes_gr)))] 
        gr_all_list[[region]] = gr_all_list[[region]][unique(from(findOverlaps(gr_all_list[[region]], wes_gr)))] 
    }
}

###################################################
################# Calc CNV corr ###################
###################################################
merged_list = c(gr_list)
if (T) {
# if (!file.exists(str_interp("${output_dir}/${case_id}__${setting}.tsv"))) {
  cl = makeCluster(16)
  registerDoParallel(cl)
  region_id_list = names(gr_list)

  cor_df = foreach(region1=region_id_list, .combine=bind_rows, .packages=c("tidyverse","GenomicRanges") ) %dopar% {
    print(region1)
    hmm_gr1 = gr_list[[region1]]
    hmm_gr_pred1 = gr_pred_list[[region1]]
    hmm_gr_all1 = gr_all_list[[region1]]
    this_cor_df = data.frame()
    for (region2 in region_id_list) {
	hmm_gr2 = gr_list[[region2]]
	hmm_gr_pred2 = gr_pred_list[[region2]]
        hmm_gr_all2 = gr_all_list[[region2]]
	jac = CNV_profile_jaccard(hmm_gr1, hmm_gr2)
	jac_pred = CNV_profile_jaccard_prob(hmm_gr_pred1, hmm_gr_pred2, mode="pred")
	jac_all = CNV_profile_jaccard_prob(hmm_gr_all1, hmm_gr_all2, mode="all")
	this_cor_df = bind_rows(this_cor_df, data.frame(id1=region1, id2=region2, jac=jac, jac_pred=jac_pred, jac_all=jac_all))
    }
    this_cor_df
}
cor_df %>% write_tsv(str_interp("${output_dir}/${case_id}__${setting}.tsv"))
}

cor_df = read_tsv(str_interp("${output_dir}/${case_id}__${setting}.tsv"))
region_id_list = unique(cor_df$id1) %>% sort
piece_id_list = gsub("__.*","",region_id_list) %>% unique %>% sort
piece_col = rainbow_hcl(length(piece_id_list)) %>% `names<-`(piece_id_list)

jac_mat = cor_df %>%  
  pivot_wider(id_cols="id2", names_from = "id1", values_from = jac) %>% 
  column_to_rownames("id2") %>% 
  as.matrix %>% 
  .[region_id_list, region_id_list]

jac_pred_mat = cor_df %>% 
  pivot_wider(id_cols="id2", names_from = "id1", values_from = jac_pred) %>% 
  column_to_rownames("id2") %>% 
  as.matrix %>%
  .[region_id_list, region_id_list]

jac_all_mat = cor_df %>% 
  pivot_wider(id_cols="id2", names_from = "id1", values_from = jac_all) %>% 
  column_to_rownames("id2") %>% 
  as.matrix %>% 
  .[region_id_list, region_id_list]

hclust_jac = hclust(as.dist(1-jac_mat), method=linkage_method)
hclust_jac_pred = hclust(as.dist(1-jac_pred_mat), method=linkage_method)
hclust_jac_all = hclust(as.dist(1-jac_all_mat), method=linkage_method)

ht=Heatmap(
  t(jac_pred_mat),
  name="Prob\nJac", 
#  col = colorRamp2(c(0, 1, 2), hcl_palette = "RdBu", rev=F), 
  col = colorRamp2(c(0, 1), hcl_palette = "Rocket", rev=T), 
  width = unit(60, "mm"), 
  height = unit(60, "mm"),
  column_title = "Predicted states",
  cluster_row_slices =F, cluster_column_slices=F,
  cluster_rows = hclust_jac, cluster_columns= hclust_jac,
  row_names_gp = gpar(fontsize = max(8, 120/nrow(jac_pred_mat))),
  column_names_gp = gpar(fontsize = max(8, 120/ncol(jac_pred_mat))),
  show_row_names = F, show_column_names =F,
  top_annotation = HeatmapAnnotation(piece=gsub("__.*","",region_id_list), show_legend=F, col=list(piece=piece_col))
)

ht2=Heatmap(
  t(jac_all_mat),
  name="Prob\nJac", 
#  col = colorRamp2(c(0, 1, 2), hcl_palette = "RdBu", rev=F), 
  col = colorRamp2(c(0, 1), hcl_palette = "Rocket", rev=T), 
  width = unit(60, "mm"), 
  height = unit(60, "mm"),
  column_title = "All states",
  cluster_row_slices =F, cluster_column_slices=F, 
  cluster_rows = hclust_jac_pred, cluster_columns= hclust_jac_pred,
  row_names_gp = gpar(fontsize = max(8, 120/nrow(jac_all_mat))),
  column_names_gp = gpar(fontsize = max(8, 120/ncol(jac_all_mat))),
  show_row_names = F, show_column_names =F,
  top_annotation = HeatmapAnnotation(piece=gsub("__.*","",region_id_list), show_legend=F, col=list(piece=piece_col))
)

ht3=Heatmap(
  t(jac_mat),
  name="Jac", 
#  col = colorRamp2(c(0, 1, 2), hcl_palette = "RdBu", rev=F), 
  col = colorRamp2(c(0, 1), hcl_palette = "Rocket", rev=T), 
  width = unit(60, "mm"), 
  height = unit(60, "mm"),
  column_title = "No Prob",
  cluster_row_slices =F, cluster_column_slices=F, 
  cluster_rows = hclust_jac_all, cluster_columns= hclust_jac_all,
  row_names_gp = gpar(fontsize = max(8, 120/nrow(jac_all_mat))),
  column_names_gp = gpar(fontsize = max(8, 120/ncol(jac_all_mat))),
  show_row_names = F, show_column_names =F,
  top_annotation = HeatmapAnnotation(piece=gsub("__.*","",region_id_list), show_legend=F, col=list(piece=piece_col))
)

ht_col=Heatmap(
  t(jac_pred_mat),
  name="Prob\nJac", 
#  col = colorRamp2(c(0, 1, 2), hcl_palette = "RdBu", rev=F), 
  col = colorRamp2(c(0, max(jac_pred_mat)), hcl_palette = "Rocket", rev=T), 
  width = unit(60, "mm"), 
  height = unit(60, "mm"),
  column_title = "Predicted states (color scale)",
  cluster_row_slices =F, cluster_column_slices=F,
  clustering_distance_rows = function(x){as.dist(1-x)}, clustering_distance_columns = function(x){as.dist(1-x)}, 
  cluster_rows = hclust_jac_pred, cluster_columns= hclust_jac_pred,
  row_names_gp = gpar(fontsize = max(8, 120/nrow(jac_pred_mat))),
  column_names_gp = gpar(fontsize = max(8, 120/ncol(jac_pred_mat))),
  show_row_names = F, show_column_names =F,
  top_annotation = HeatmapAnnotation(piece=gsub("__.*","",region_id_list), show_legend=F, col=list(piece=piece_col))
)

ht2_col=Heatmap(
  t(jac_all_mat),
  name="Prob\nJac", 
#  col = colorRamp2(c(0, 1, 2), hcl_palette = "RdBu", rev=F), 
  col = colorRamp2(c(0, max(jac_all_mat)), hcl_palette = "Rocket", rev=T), 
  width = unit(60, "mm"), 
  height = unit(60, "mm"),
  column_title = "All states (color scale)",
  cluster_row_slices =F, cluster_column_slices=F, 
  clustering_distance_rows = function(x){as.dist(1-x)}, clustering_distance_columns = function(x){as.dist(1-x)}, 
  cluster_rows = hclust_jac_all, cluster_columns= hclust_jac_all,
  row_names_gp = gpar(fontsize = max(8, 120/nrow(jac_all_mat))),
  column_names_gp = gpar(fontsize = max(8, 120/ncol(jac_all_mat))),
  show_row_names = F, show_column_names =F,
  top_annotation = HeatmapAnnotation(piece=gsub("__.*","",region_id_list), show_legend=F, col=list(piece=piece_col))
)

p_hist = ggplot(data.frame(jac_pred=as.vector(jac_pred_mat)), aes(x=jac_pred)) +
	geom_histogram() +
	theme_classic() +
	labs(title="Jac (predicted state)")

p_list = list("jac_pred_ht"=grid.grabExpr(draw(ht)), 
	      "jac_all_ht"=grid.grabExpr(draw(ht2)), 
	      "jac_ht"=grid.grabExpr(draw(ht3)), 
	      "jac_pred_ht_col"=grid.grabExpr(draw(ht_col)), 
	      "jac_all_ht_col"=grid.grabExpr(draw(ht2_col)),
	      "jac_pred_hist"=p_hist)

pdf(str_interp("${output_dir}/${case_id}__${setting}.pdf"), w=12, h=8)
wrap_plots(p_list, nrow=2)
dev.off()

st_tumor_list = list()
for (sample_id in sample_id_list) {
    st = read_st_obj(sample_id)
    subclone_tbl = read_subclone_tbl(sample_id)
    st = AddMetaData(st, subclone_tbl[Cells(st),])
    st_tumor = subset(st, Filtered_tumor_regions!="0")
    st_tumor_list[[sample_id]] = st_tumor
}

p_st_list = list()
for (sample_id in sample_id_list) {
    st_tumor = st_tumor_list[[sample_id]]
    p = SpatialPlot(st_tumor, group.by="Filtered_tumor_regions", stroke=0, label=T, crop=F, label.box=F, label.size=3) + 
	    scale_fill_discrete(na.value=NA) + labs(title=sample_id) + NoLegend()
    p_st_list[[sample_id]] = p
}

distmat_selected = jac_pred_mat
hclust_selected = hclust(as.dist(1-distmat_selected), method=linkage_method)
opt_k = maptree::kgs(hclust_selected, as.dist(1-distmat_selected), alpha=2, maxclus=MAX_K)
opt_k[which(opt_k == min(opt_k))]

for (K in c(2:MAX_K)) {
    all_clust_id = hclust_selected %>% cutree(k=K) 
    for (sample_id in sample_id_list) {
	clust_col = rainbow_hcl(K) %>% `names<-`(seq(K))

	st_tumor = st_tumor_list[[sample_id]] 
	sample_clust_id = all_clust_id[grep(sample_id, names(all_clust_id), value=T)]
        names(sample_clust_id) = gsub(".*__","",names(sample_clust_id))
	st_tumor$clust_id = sample_clust_id[as.character(st_tumor$Filtered_tumor_regions)] %>% as.character
	p = SpatialPlot(st_tumor, group.by="clust_id", crop=F, image.alpha=0) + 
	    scale_fill_manual(values=clust_col, na.value=NA) + labs(title=str_interp("K=${K}")) + NoLegend()
        p_st_list[[paste0(sample_id,"__", K)]] = p
    }
}

png(str_interp("${output_dir}/Clone_${case_id}.png"), h=2*length(st_tumor_list), w=2*MAX_K, units="in", res=150)
wrap_plots(p_st_list, ncol=MAX_K, byrow=F)
dev.off()
############ End of Calc CNV corr ###################
