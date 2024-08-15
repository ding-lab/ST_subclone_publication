# This script find spot from adjacent slice that is closest to the spot in the current slice
library(tidyverse)
library(patchwork)
library(googlesheets4)
library(tictoc)
library(future)
library(furrr)

plan(multisession, workers = 10)

# Out
project_path = ""
out_path =str_glue("{project_path}/2_PASTE2_anaysis/2_get_matching_spot/out/")
dir.create(out_path, recursive = T, showWarnings = F)

# Load tracking sheet 
gs4_deauth()
google_sheet_url = ""
sheet_multi = read_sheet(google_sheet_url, sheet = 'Sample_multisections')

# Load subclone
subclone_df = read_tsv(str_glue('{project_path}/0_data/table/OCT_genetic_clones_2023-09-18.tsv'))

# List all case/piece to plot
paste2_result_root = str_glue('{project_path}/1_run_PASTE2/out')
all_piece_names = paste2_result_root %>% {intersect(list.files(.), list.dirs(., full.names= F))}

for(piece_use in all_piece_names){
    ## ------------------ 1. Combine all alignment result and meta.data ------------------ ###
    # 1. Load PASTE2 coordinates
    message(piece_use)
    coord_path = str_glue("{paste2_result_root}/{piece_use}")

    sample_names = list.files(coord_path, pattern = 'csv') %>% str_remove('.csv') %>% gtools::mixedsort(.)

    all_coord = map(sample_names, function(sample_nm){
        read_csv(str_glue("{coord_path}/{sample_nm}.csv")) %>% 
            mutate(
                PASTEpieceName = piece_use,
                PASTEFileName = sample_nm
                ) 
    }) %>% bind_rows()

    # 2. Load morph regions
    sheet_piece_use = sheet_multi %>% filter(PASTEpieceName== piece_use)
    morph_df = pmap(sheet_piece_use, ~with(list(...), {
            read_tsv(MorphOutput) %>% 
            dplyr::rename(barcode = '...1') %>%
            # Clean columns
            select(barcode, `Filtered tumor regions`) %>%
            dplyr::rename(Filtered_tumor_regions= `Filtered tumor regions`) %>%
            mutate(
                PASTEFileName = all_of(PASTEFileName),
                PASTEpieceName = all_of(PASTEpieceName),
                LibraryName = all_of(LibraryName)
                )
        })) %>% 
        bind_rows() 
    # Add morph regions to coord
    all_coord = left_join(all_coord, morph_df, by = c('PASTEpieceName', 'PASTEFileName', 'barcode'))

    # 3. Add genetic subclone and libraryName
    sheet_multi_use = sheet_multi %>% filter(PASTEpieceName == piece_use) %>% select(PASTEpieceName, PASTEFileName, LibraryName)
    all_coord_subclone = left_join(all_coord, subclone_df, c('LibraryName' = 'sample_id','Filtered_tumor_regions'))
    
    # ------------------------------------------ ###
    # Make connection between each spot to its adjacent slice
    # 4. Add 
    # Generate pariwise plot
    sheet_piece_use = sheet_multi %>% filter(PASTEpieceName == piece_use)
    piece_name = sheet_piece_use$PASTEpieceName %>% unique
    # number_range = str_remove(sample_names, piece_name) %>% as.numeric %>% sort
    n_slice = length(sample_names)
    plot_group = data.frame(
        section1 = sample_names[1:(n_slice-1)],
        section2 = sample_names[2:n_slice]
    )

    # method 1b, calculate min distance. take longer
    # Rduce search space to +- 100
    tic('Outer loop')
    result_list = pmap(plot_group, function(section1, section2){
        tic('Inner loop')
        message(section1, ' and ',section2)
        section1_df = all_coord_subclone %>% filter(PASTEFileName == section1) %>% setNames(str_c(names(.), '_1'))
        section2_df = all_coord_subclone %>% filter(PASTEFileName == section2) %>% setNames(str_c(names(.), '_2'))

        section1_2_df = section1_df %>% 
            future_pmap(~with(list(...), {
            section2_df %>% 
                # Add section 1 info - spot level
                mutate(
                    #x_1 = all_of(x_1), y_1 = all_of(y_1),
                    barcode_1 = all_of(barcode_1),
                    Filtered_tumor_regions_1 = all_of(Filtered_tumor_regions_1), 
                    genetic_clone_1 = all_of(genetic_clone_1),
                    PASTEpieceName_1 = all_of(PASTEpieceName_1),
                    PASTEFileName_1 = all_of(PASTEFileName_1),
                    LibraryName_1 = all_of(LibraryName_1)
                ) %>%
                # reduce search space to only +- 100
                filter(between(x_2, x_1 - 100, x_1 + 100), between(y_2, y_1 - 100, y_1 + 100)) %>%
                mutate(spot_dist = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)) %>%
                slice(which.min(spot_dist)) %>% 
                mutate(x_1 = all_of(x_1), y_1 = all_of(y_1))
        })) %>% bind_rows() %>%
            # section level info
            mutate(section1 = section1, section2 = section2) 
        # Reorder to put all _1 to the front
        section1_2_df = section1_2_df %>% select(section1, section2, spot_dist, ends_with('_1'), ends_with('_2'))
        toc(); print(Sys.time())
        return(section1_2_df)
    })
    toc(); print(Sys.time())

    # Save result
    #dir.create(str_glue("{out_path}"), recursive = T, showWarnings = F)
    result_list %>% bind_rows() %>% 
        write_csv(str_glue("{out_path}/{piece_name}_min_dist.csv"))
    
}



