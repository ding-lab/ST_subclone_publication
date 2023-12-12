# Run RCTD with matching snRNA
# Jingxian Clara Liu
# 01/17/23

suppressMessages(library(Seurat))
# suppressMessages(library(Signac))
suppressMessages(library(tidyverse))
suppressMessages(library(sctransform))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))
suppressMessages(library(future))
suppressMessages(library(RCTD))
suppressMessages(library(Matrix))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scatterpie))
suppressMessages(library(paletteer))
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))

RhpcBLASctl::blas_set_num_threads(20)
# plan("multiprocess", workers = 10)
options(future.globals.maxSize = 250 * 1024 ^ 3)
options(future.fork.multithreading.enable = FALSE)
plan()
set.seed(1234)

source("/diskmnt/Projects/Users/cliu/pancan_ST/Overview/helper_global.R")
selected_tbl = read_tracking_sheet()

### Some constants --------------------------------------------------
# Minimal number of cells per cell type in the snRNA obj to be counted in RCTD reference
MIN_CELL_COUNT_RCTD_REF = 50

### Set sample id to be processed -----------------------------------
sample_id <- "HT260C1-Th1K1Fc2U1Z1Bs1"
sample_id <- "HT308B1-S1H1Fc2U1Z1Bs1"
output_dir <- paste0("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/13-RCTD/13_1-RCTD_output/current/Sample_ref/")
celltype_colname = "cell_type_final"

manual_celltype_all = read_tsv("/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/28-snRNA/28_0-snRNA_obj/sn_list_metadata_manual.tsv")
unique_snRNA_ids = selected_tbl$Paired_snRNA_id %>% unique %>% setdiff(c("No_snRNA",NA))
manual_celltype_all %>% 
	filter(Paired_snRNA_id %in% unique_snRNA_ids) %>% 
	distinct(Paired_snRNA_id, cell_type_final) %>%
	pull(cell_type_final) %>% table %>% sort %>% rev

### Parameters ------------------------------------------------------
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

### RCTD with matching snRNA ----------------------------------------
# if (sn_obj_path=="") {sn_obj_path = selected_tbl[sample_id, "Paired_snRNA_obj"]}
# if (st_obj_path=="") {st_obj_path = selected_tbl[sample_id, "SeuratObject"]}

sn_obj_path = selected_tbl[sample_id, "Paired_snRNA_obj"]
snRNA_id = selected_tbl[sample_id, "Paired_snRNA_id"]
snRNA_meta = selected_tbl[sample_id, "Paired_snRNA_metadata"] 

st <- read_st_obj(sample_id)
sn = readRDS(sn_obj_path)

sample_output_dir = str_interp("${output_dir}/${sample_id}/")
if (!dir.exists(sample_output_dir)) {dir.create(sample_output_dir)}

snRNA_meta_tbl = read_tsv(snRNA_meta) %>% 
	column_to_rownames("barcode")
sn = AddMetaData(sn, snRNA_meta_tbl[Cells(sn), ])

# Remove doublets/unknowns
cells.use = Cells(sn)[!(sn@meta.data[,celltype_colname] %in% 
			c("Doublet","DOUBLET","doublet","Remove_Unknown","Remove_badSample"))]
sn <- subset(sn, cells = cells.use)

# Change T_reg to T
sn@meta.data[,celltype_colname] = gsub("_reg","",sn@meta.data[,celltype_colname])

# Load reference snRNA data
celltype_count = table(sn@meta.data[,celltype_colname])
cell_types_big <- celltype_count[celltype_count > MIN_CELL_COUNT_RCTD_REF] %>% names()
cell_types_big

cells.use = Cells(sn)[sn@meta.data[,celltype_colname] %in% cell_types_big ]
sn <- subset(sn, cells = cells.use)

counts <- sn[["RNA"]]@counts
cell_types <- factor(sn@meta.data[,celltype_colname], cell_types_big)
names(cell_types) = Cells(sn)
sn.ref <- Reference(counts, cell_types)

print("\n\n[sn Refernce Cell type]\n")
table(sn@meta.data[,celltype_colname])
print("\n[sn Refernce Cell type]\n")

# Load query spatial transcriptome data
counts <- st[["Spatial"]]@counts
coords <- as.data.frame(GetTissueCoordinates(st))
coords <- coords[,c(2,1)]
coords[,2] <- -coords[,2]
st.query <- SpatialRNA(coords, counts)

# Plot nUMI
st.barcodes <- colnames(st.query@counts)
pdf(paste0(sample_output_dir, "/RCTD_",sample_id,"_nUMI.pdf"))
plot_puck_continuous(st.query, st.barcodes, st.query@nUMI, ylimit = c(0,round(quantile(st.query@nUMI,0.9))), 
                     title ='plot of nUMI')
dev.off()

# run RCTD
myRCTD <- create.RCTD(st.query, sn.ref, max_cores = 1)
# MODE = ("doublet","full", "multi")
# mode = "doublet"
# myRCTD <- run.RCTD(myRCTD, doublet_mode = mode)
# saveRDS(myRCTD, paste0(paste0(sample_output_dir, '/RCTD_',mode,"_",sample_id,'.rds')))
mode = "multi"
myRCTD <- run.RCTD(myRCTD, doublet_mode = mode)
saveRDS(myRCTD, paste0(paste0(sample_output_dir, '/RCTD_',mode,"_",sample_id,'.rds')))
# mode = "full"
# myRCTD <- run.RCTD(myRCTD, doublet_mode = mode)
# saveRDS(myRCTD, paste0(paste0(sample_output_dir, '/RCTD_',mode,"_",sample_id,'.rds')))

## Default plotting for full mode:
if (mode=="full") {
  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA 

  # Plots the confident weights for each cell type as in full_mode ('cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, spatialRNA, sample_output_dir, norm_weights) 
  # Plots all weights for each cell type as in full_mode ('cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, sample_output_dir, norm_weights) 
  # Plots the number of confident pixels of each cell type in 'full_mode' ('cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, sample_output_dir, norm_weights, spatialRNA)
}
