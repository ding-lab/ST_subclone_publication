# conda activate seurat4.3
library(optparse)
library(tidyverse)
library(patchwork)
library(googlesheets4)
library(tictoc)
library(future)
library(furrr)
library(igraph)
library(ggraph)
library(Seurat)

# 1. Overall level: number of component, aka number of the 3D volume 
# 2. 3D volume level: number of cycle, handle
# 3. 3D volume level: maximum degree of each volume (espeically for sample with only 2 sections since there will be no loops)

# color
project_folder = ""
analysis_folder = str_glue('{project_folder}/2_PASTE2_analysis/5_calculate_metrics')
source(str_glue('{project_folder}/0_palette/color_STsubclone.r'))

# parameters
out_path = str_glue("{analysis_folder}/out")
dir.create(out_path, recursive = T, showWarnings = F)


######## 0. Loading data and info ########################################################
# Load tracking sheet 
gs4_deauth()
google_sheet_url = ""
sheet_multi = read_sheet(google_sheet_url, sheet = 'Sample_multisections') %>% 
    filter(!is.na(PASTEFileName)) 


# Load inidvidual PASTE2 coordinates
paste2_df_list = pmap(sheet_multi, ~with(list(...), {
    message(LibraryName,' ' , piece_id, ' ' , PASTEFileName)
    read_csv(PASTE2XinHaoCoordinate) %>% 
        mutate(
            PASTEpieceName = PASTEpieceName,
            PASTEFileName = PASTEFileName,
            LibraryName = LibraryName
        )
})) %>% setNames(sheet_multi$LibraryName) 

# Load objects
obj_list = pmap(sheet_multi, ~with(list(...), {
    message("Loading :", LibraryName, "Object")
    read_rds(SeuratObject)
})) %>% setNames(sheet_multi$LibraryName)

# # Load region connection 
three_d_cluster_root = str_glue("{analysis_folder}/4_multi_section_regions_Ver20231031_N_min_3")
all_pieces = list.files(str_glue('{three_d_cluster_root}/out'))
three_d_cluster_list = map(all_pieces, ~{
    read_tsv(str_glue('{three_d_cluster_root}/out/{.x}/4a_3D_cluster.tsv'))
}) %>% setNames(all_pieces)

## ----- ANALYSIS  ------------------------------------------------------------
# 1. Create graph 
message("Using undirected graph")
graph_list = map(three_d_cluster_list, function(three_d_cluster_df_use){
    # Split by each volumn
    g_list = three_d_cluster_df_use %>% split(., f = .$volumn_cluster1) %>% 
        map(function(df, volumn_id){
            edges_df = df %>% select(region1, region2) %>% setNames(c('From','To'))
            nodes_df = data.frame(name = union(df$region1, df$region2))
            
            g <- graph_from_data_frame(edges_df, directed=FALSE, vertices=nodes_df)
            return(g)
        })
}) %>% unlist(recursive = F) %>% as.list() 

## FUNCTIONS ------------------------------------------------------------

#### NEW 2023/10/31
count_loop = function(graph_obj){
    # n loop = n-edges - n-nodes + 1
    length(E(graph_obj)) - length(V(graph_obj)) + 1
}

plot_degree_distribution = function(graph_obj, title = 'Degree distribution'){
    data.frame(
        degree_frequency = degree_distribution(graph_obj)
    ) %>% mutate(degree = 0:(nrow(.)-1)) %>% 
    ggplot(aes(x = degree, y = degree_frequency, fill = degree_frequency)) + 
        geom_bar(stat = 'identity') + 
        cowplot::theme_cowplot() +
        labs(title = title, x = 'Degree', y = 'Frequency') +
        scale_fill_viridis_c()
}

mlutiple_degree_distribution_df = function(graph_list){
    imap(graph_list, ~{
        data.frame(
            degree_frequency = degree_distribution(.x)
        ) %>% mutate(degree = 0:(nrow(.)-1)) %>% 
        mutate(volume_id = .y)
    }) %>% bind_rows()
}
plot_mlutiple_degree_distribution = function(graph_list){
    mlutiple_degree_distribution_df(graph_list) %>% 
    ggplot(aes(x = degree, y = degree_frequency, fill = degree)) + 
        geom_bar(stat = 'identity') + 
        facet_wrap(~volume_id) + 
        theme_bw() +
        labs(x = 'Degree', y = 'Frequency')+
        scale_fill_viridis_c()
}

## ----- ANALYSIS  ------------------------------------------------------------
# Calculate metrics
### Metric: LOOPS
analysis_loop_out = str_glue('{out_path}/0_loop_per_volume/');dir.create(analysis_loop_out, recursive = T, showWarnings = F)
loop_per_volume_df = map(graph_list, count_loop) %>% unlist %>% as.data.frame %>% setNames('n_loop') %>% rownames_to_column('volume_id') %>% 
    separate(volume_id, into = c('piece', 'volume_id'), sep = '\\.', remove = F) 
write_tsv(loop_per_volume_df, str_glue('{analysis_loop_out}/0_loop_per_volume.tsv'))
# none zero entry
loop_per_volume_df %>% filter(n_loop >0)

### Metric: DEGREE
# A. Plot all degree distribution
analysis_a_out = str_glue('{out_path}/A_degree_distribution_all/');dir.create(analysis_a_out, recursive = T, showWarnings = F)
for(piece_use in all_pieces){
    graph_list_use = graph_list %>% .[str_subset(names(.), piece_use)]
    if(length(graph_list_use) == 0) next()
    distribution_df = mlutiple_degree_distribution_df(graph_list_use) 
    write_tsv(distribution_df, str_glue('{analysis_a_out}/A_degree_distribution_{piece_use}.tsv'))
    pdf(str_glue('{analysis_a_out}/0_degree_distribution_{piece_use}.pdf'))
    plot_mlutiple_degree_distribution(graph_list_use) %>% print
    dev.off()
}

# 2. Plot metrics per sample
# B. Plot maximum degree distribuiotn
analysis_b_out = str_glue('{out_path}/B_degree_distribution_max/');dir.create(analysis_b_out, recursive = T, showWarnings = F)
for(piece_use in all_pieces){
    graph_list_use = graph_list %>% .[str_subset(names(.), piece_use)]
    if(length(graph_list_use) == 0) next()
    distribution_df = mlutiple_degree_distribution_df(graph_list_use) %>%
        dplyr::rename(degree_frequency_per_volume = degree_frequency) 
    # keep max degree per each volume
    distribution_max_df = distribution_df %>% group_by(volume_id) %>% slice_max(degree) %>% 
        dplyr::rename(max_degree = degree)
    write_tsv(distribution_max_df, str_glue('{analysis_b_out}/B_degree_distribution_max_{piece_use}.tsv'))
    # plot
    pdf(str_glue('{analysis_b_out}/B_degree_distribution_max_{piece_use}.pdf'))
    p_dist_max = distribution_max_df %>%
        ggplot(aes(x = volume_id, y = max_degree, fill = max_degree)) + 
        geom_bar(stat = 'identity') + 
        theme_bw() +
        labs(x = 'volume ID', y = 'Maximum Degree per volume')#+
        #scale_fill_viridis_d()
    print(p_dist_max)
    dev.off()
}

# C. Plot maxnium degree distribution of all samples at onces!
analysis_c_out = str_glue('{out_path}/C_degree_distribution_max_all/');dir.create(analysis_c_out, recursive = T, showWarnings = F)
distribution_max_allsamples_df = map(all_pieces, function(piece_use){
    graph_list_use = graph_list %>% .[str_subset(names(.), piece_use)]
    if(length(graph_list_use) == 0) return(NULL)
    distribution_df = mlutiple_degree_distribution_df(graph_list_use) %>%
        dplyr::rename(degree_frequency_per_volume = degree_frequency) 
    # keep max degree per each volume
    distribution_max_df = distribution_df %>% group_by(volume_id) %>% slice_max(degree) %>% 
        dplyr::rename(max_degree = degree)
    return(distribution_max_df)
}) %>% bind_rows() %>% 
separate(volume_id, into = c('piece', 'volume'), sep = '\\.', remove = F) 

# c2. Add in meta data : cancer type, metastasis, n sections
sheet_multi_clean = sheet_multi %>% distinct(PASTEpieceName, N_sections, cancer_type, sample_type_tumor) %>% dplyr::rename(piece = PASTEpieceName)
distribution_max_allsamples_df = left_join(distribution_max_allsamples_df, sheet_multi_clean, by = 'piece')

# c3. Add n components
distribution_max_allsamples_df = distribution_max_allsamples_df %>% 
    group_by(piece) %>% 
    mutate(n_componenet = n()) %>% 
    # sort sample by n_component per cancer type
    arrange(cancer_type, n_componenet) 
sample_order = distribution_max_allsamples_df$piece %>% unique %>% rev
write_tsv(distribution_max_allsamples_df, str_glue('{analysis_c_out}/C3_degree_distribution_max_all.tsv'))

p_max_dist_all = distribution_max_allsamples_df %>% ggplot(aes(x = factor(piece, levels = sample_order), y = max_degree, fill = cancer_type)) + 
    facet_grid(.~cancer_type, space = 'free', scales='free') + 
    geom_violin(color = 'transparent', scale = 'width', alpha = 0.6) + 
    geom_boxplot(width = 0.1, alpha = 0.3, outlier.shape = NA) + 
    geom_jitter(aes(fill = cancer_type), color = 'gray20', height = 0, size = 3, seed = 123, alpha = 0.7, shape = 21) + 
    labs(title = 'Maximum degree distribution per sample', x = 'Sample', y = 'Maximum degree') +
    theme_bw() + NoLegend() + 
    scale_fill_manual(values = color_cancer_all) 

pdf(str_glue('{analysis_c_out}/C3_degree_distribution_max_all.pdf'),w = 10, h = 5)
p_max_dist_all %>% print
dev.off()


# c4. Plot n component per cancer type
p_n_component_cancer_type = distribution_max_allsamples_df %>% 
    distinct(piece, n_componenet, cancer_type) %>% 
    ggplot(aes(x = factor(piece, levels = sample_order), y = n_componenet, fill = cancer_type)) + 
    facet_grid(.~cancer_type, space = 'free', scales='free') + 
    geom_bar(stat = 'identity') + 
    # Add number, center
    geom_text(aes(label = n_componenet, y = n_componenet/2), vjust = -0.5, size = 3) +
    labs(title = 'Number of components (tumor volume) per sample', x = 'Sample piece', y = 'Number of components (tumor volume)') +
    theme_bw() + NoLegend() + 
    scale_fill_manual(values = color_cancer_all)

pdf(str_glue('{analysis_c_out}/C4_n_component_per_sample.pdf'),w = 10, h = 4)
p_n_component_cancer_type %>% print
dev.off()


