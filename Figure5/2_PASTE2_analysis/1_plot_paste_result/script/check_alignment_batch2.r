library(tidyverse)
library(patchwork)

# Out
project_path = ""
out_path =str_glue("{project_path}/2_PASTE2_anaysis/1_plot_paste_result/out/")

dir.create(out_path, recursive = T, showWarnings = F)

# List all case/piece to plot
paste_result_path = str_glue("{project_path}/1_run_PASTE2/out")
all_piece_names = paste_result_path %>% {intersect(list.files(.), list.dirs(., full.names= F))}

for(piece_use in all_piece_names){
    # Coord path
    coord_path = str_glue("{paste_result_path}/{piece_use}")

    sample_names = list.files(coord_path, pattern = 'csv') %>% str_remove('.csv')

    all_coord = map(sample_names, function(sample_nm){
        read_csv(str_glue("{coord_path}/{sample_nm}.csv")) %>% 
            mutate(sample_name = sample_nm) 
    }) %>% bind_rows()


    # Generate pariwise plot
    piece_name = str_remove(sample_names, '[0-9]$') %>% unique
    number_range = str_remove(sample_names, piece_name) %>% as.numeric %>% sort
    n_slice = length(sample_names)
    plot_group = data.frame(
        section1 = str_c(piece_name, seq(min(number_range), max(number_range-1))),
        section2 = str_c(piece_name, seq(min(number_range) + 1, max(number_range)))
    )


    p_pariwise = pmap(plot_group, function(section1, section2){
        all_coord %>% filter(sample_name %in% c(section1, section2)) %>% 
            ggplot(aes(x = x ,y =y, color = sample_name)) + 
            geom_point(size = 0.7) + 
            labs(title = str_glue("{section1} vs {section2}")) +
            theme_bw() + theme(aspect.ratio = 1) 
    })

    p_all = ggplot(all_coord, aes(x = x ,y =y, color = sample_name)) + 
        geom_point(size = 0.7) +
        theme_bw() + theme(aspect.ratio = 1) + 
        labs(title = "All sections")

    # put in seperate folder
    outpath_piece = str_glue("{out_path}/{piece_use}/")
    dir.create(outpath_piece, recursive = T, showWarnings = F)
    pdf(str_glue("{outpath_piece}/{piece_name}_all_and_pairwise.pdf"), width = 10, height = 10)
    wrap_plots(c(list(p_all), p_pariwise), ncol = 2) %>% print()
    dev.off()

}


### ------------------ 2. Check alignment result ------------------ ###
library(googlesheets4)
# Pareters 
piece_name
# Part 2 uses seurat clustesr to check alignment result
gs4_deauth()
googlt_sheet_url=""
sheet_use = read_sheet(googlt_sheet_url, sheet = 'Sample_multisections') %>% 
    filter(PASTEBatch == 'Batch2')

 # Load seurat objects
obj_list = pmap(sheet_use, function(SeuratObject,...){
        readRDS(SeuratObject)
    }) %>% setNames(sheet_use$PASTEName)

# Test loading sample
for(piece_use in all_piece_names){
    message(piece_use)
    # Load object
    sheet_piece_use = sheet_use %>% filter(PASTEpieceName == piece_use)
    obj_list_use = pmap(sheet_piece_use, function(SeuratObject,...){
        readRDS(SeuratObject)
    }) %>% setNames(sheet_piece_use$PASTEFileName)

    # Coord path
    coord_path = str_glue("{paste_result_path}/{piece_use}")

    sample_names = list.files(coord_path, pattern = 'csv') %>% str_remove('.csv')

    all_coord = map(sample_names, function(sample_nm){
        read_csv(str_glue("{coord_path}/{sample_nm}.csv")) %>% 
            mutate(sample_name = sample_nm) 
    }) %>% bind_rows()
    
    # Add seurat clusters to table
    all_clusters = obj_list_use %>% imap(function(obj, sample_nm){
        obj@meta.data %>% rownames_to_column('barcode') %>% 
            select(barcode, seurat_clusters) %>%
            mutate(sample_name = sample_nm)
    }) %>% bind_rows()

    all_coord_clusters = left_join(all_coord, all_clusters, by = c('sample_name', 'barcode'))
   
    # PLOT
    # Generate pariwise plot
    piece_name = sheet_piece_use$PASTEpieceName %>% unique
    # number_range = str_remove(sample_names, piece_name) %>% as.numeric %>% sort
    n_slice = length(sample_names)
    plot_group = data.frame(
        section1 = sample_names[1:(n_slice-1)],
        section2 = sample_names[2:n_slice]
    )

    p_pariwise = pmap(plot_group, function(section1, section2){
        all_coord_clusters %>% 
            filter(sample_name %in% c(section1, section2)) %>% 
            ggplot(aes(x = x ,y =y, color = seurat_clusters)) + 
            geom_point(size = 0.7, alpha = 0.4) +
            labs(title = str_glue("{section1} vs {section2}")) +
            #scale_color_manual(values = c('black', 'red')) +
            theme_bw() + theme(aspect.ratio = 1) 
    })
    p_individuals = map(sample_names, function(slice_use){
        all_coord_clusters %>% filter(sample_name %in% c(slice_use)) %>% 
            ggplot(aes(x = x ,y =y, color = seurat_clusters)) + 
            geom_point(size = 0.7) + 
            labs(title = str_glue("{slice_use}")) +
            theme_bw() + theme(aspect.ratio = 1) 
    })

    p_all = ggplot(all_coord_clusters, aes(x = x ,y =y, color = seurat_clusters)) + 
        geom_point(size = 0.7, alpha = 0.4) + 
        theme_bw() + theme(aspect.ratio = 1) + 
        labs(title = "All sections")

    # put in seperate folder
    outpath_piece = str_glue("{out_path}/{piece_use}/")
    dir.create(outpath_piece, recursive = T, showWarnings = F)
    pdf(str_glue("{outpath_piece}/{piece_name}_all_and_pairwise_seurat_clusters.pdf"), width = 10, height = 10)
    wrap_plots(c(list(p_all), p_pariwise), ncol = 2) %>% print()
    wrap_plots(p_individuals, ncol = 2) %>% print()
    dev.off()

}
