# SEARCH CRITERIA for the definition
library(tidyverse)
library(patchwork)
library(googlesheets4)
library(tictoc)
library(future)
library(furrr)
library(igraph)

# Functions
# Create Dummy Sankey Link
source('PATH/TO/ANALYSIS/FOLDER/4_multi_section_regions/script/src/function_sankey_make_dummy_link.r')
# Color 
source('PATH/TO/ANALYSIS/FOLDER/4_multi_section_regions/script/src/function_sankey_colors.r')


## optparse list to input necessory info
library(optparse)
option_list = list(
    make_option(c("-a", "--analysis_folder"), type="character", default=NULL, 
                help="analysis_folder"),
    make_option(c("-p", "--piece_use"), type="character", default=NULL, 
                help="piece_use"),
    make_option(c("-m", "--matching_out_path"), type="character", default=NULL, 
                help="matching_out_path"),
    # Add n_spot threshold
    make_option(c("-n", "--n_spot_minimum"), type="numeric", default=NULL, 
                help="n_spot_minimum")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# parameters
piece_use = opt$piece_use
n_spot_minimum = opt$n_spot_minimum 

analysis_folder = opt$analysis_folder
out_path =str_glue("{analysis_folder}/out/{piece_use}")
dir.create(out_path, recursive = T, showWarnings = F)


######## 0. Loading data and info ########################################################
# Load tracking sheet 
gs4_deauth()
sample_info_url = ""
sheet_multi = read_sheet(sample_info_url, sheet = 'Sample_multisections') %>% 
    filter(!is.na(PASTEFileName)) %>% 
    filter(PASTEpieceName == piece_use)


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

# Load matching spots
matching_out_path=opt$matching_out_path
if(!file.exists(str_glue('{matching_out_path}/{piece_use}_min_dist.csv'))){
    message(str_glue('{matching_out_path}/{piece_use}_matching_spot.tsv not exist'))
    message("Quit this run.")
    quit()
}
matching_df = read_csv(str_glue('{matching_out_path}/{piece_use}_min_dist.csv'))

# Load objects
obj_list = pmap(sheet_multi, ~with(list(...), {
    message("Loading :", LibraryName, "Object")
    read_rds(SeuratObject)
})) %>% setNames(sheet_multi$LibraryName)

# # 1. create networks of spots on the sample slide
##############################################################################################
##############################################################################################
# 1. Count connections between slices
# 1a. Get adjacent spots from the same sections
# NEW: 2023/10/31 Filter out connection less than n_spot_minimum spots
count_connection_df = matching_df %>% 
    # filter out NonTumor spots
    filter(
        Filtered_tumor_regions_1 != '0',
        Filtered_tumor_regions_2 != '0'
        ) %>% 
    count(section1, section2, Filtered_tumor_regions_1, Filtered_tumor_regions_2) %>%
    # v BUG: THIS SOMEHOW DOESN"T WORK. NOT SURE WHY. 
    filter(n > n_spot_minimum) %>%  # Filter out connection less than n spots
    group_by(Filtered_tumor_regions_1) %>%
    mutate(percent = n/sum(n)) %>%
    mutate(
    region1 = str_c(section1,"-",Filtered_tumor_regions_1),
    region2 = str_c(section2,"-",Filtered_tumor_regions_2)
    )

write_tsv(count_connection_df, str_glue('{out_path}/1_count_connection.tsv'))

# 2. Make sankey
### Functions ----------------------------------------------------------------------------------------
# Set color in sankey Javascript format
library(colorspace)
define_colors = function(ident, colors){
    ident_string = ident %>%  paste(collapse = "','")
    colors_string = colors %>%  paste(collapse = "','")
    my_color <- str_glue("d3.scaleOrdinal() .domain(['{ident_string}']) .range(['{colors_string}'])")
    return(my_color)
}
my_color = define_colors(c('Th1H3U1','Th1H3U2'), c('#1f77b4','#ff7f0e'))


# Library
library(networkD3)
library(dplyr)
library(gtools)

links = count_connection_df  %>% mutate(value = n)

# Works but need iteratively add dummy levels for each level
missing_2nd_level = setdiff(
    links %>% filter(section1 == 'Th1H3U2') %>% pull(region1),
    links %>% filter(section1 == 'Th1H3U1') %>% pull(region2)
    )
links = links %>% ungroup %>% add_row(region1 = 'Dummy_level1', region2 = missing_2nd_level, value = 0.1)
#### --------------
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$region1), 
  as.character(links$region2)) %>% unique()
)
 
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$region1, nodes$name)-1 
links$IDtarget <- match(links$region2, nodes$name)-1


# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              #iterations = 0 ,
              colourScale = my_color,
              sinksRight=FALSE)
p

# Save the result as html
# https://stackoverflow.com/questions/56086341/how-to-save-simplenetwork-output-from-networkd3-in-pdf-jpeg-tiff-format  
require(htmlwidgets)
saveWidget(p, file=str_glue("{out_path}/2_sankey.html"))

require(webshot) # need conda install phantomjs
webshot(url = str_glue("{out_path}/2_sankey.html"), 
    file= str_glue("{out_path}/2_sankey.pdf"),
    vwidth = 700,
    vheight = 400,
)



###### 3. Make spatial plot to check distribution ----------------------------------------------
# 3. color by micoreiongion
spatial_map_df = bind_rows(
    matching_df[,c('x_1','y_1','Filtered_tumor_regions_1','section1')] %>% rename(x=x_1, y=y_1, Filtered_tumor_regions=Filtered_tumor_regions_1, section=section1),
    matching_df[,c('x_2','y_2','Filtered_tumor_regions_2','section2')] %>% rename(x=x_2, y=y_2, Filtered_tumor_regions=Filtered_tumor_regions_2, section=section2)
    ) %>% 
    filter(
        Filtered_tumor_regions != '0'
        ) %>%
    mutate(section = factor(section, levels = mixedsort(unique(section))))
write_tsv(spatial_map_df, str_glue('{out_path}/3_Microregion_spatial_map.tsv'))
spatial_label_df = spatial_map_df %>% group_by(section,Filtered_tumor_regions) %>% summarize(across(where(is.numeric), mean))
p_microregion = ggplot(spatial_map_df, aes(x=x,y=-y, color=as.character(Filtered_tumor_regions))) +
    facet_wrap(~section) +
    geom_point(size = 0.5) + 
    geom_text(data = spatial_label_df, aes(label = Filtered_tumor_regions), color = 'gray20') + 
    cowplot::theme_cowplot() +
    theme(aspect.ratio = 1) +
    theme(
        legend.position = "bottom"
        ) + 
    labs(title = 'Microregion')
ggsave(str_glue("{out_path}/3_Microregion_spatial_map.pdf"), p_microregion, width = 10, height = 5, units = "in", dpi = 300)
    
    
############################################################################
# 4. 3D cluster
# 4a Assign 3D cluster
# CRITERIA:
# Rule of connection: if > n_spot_minimum spots are connected 
# OR > 20% are connected, they are in the same cluster
connected_df = count_connection_df %>%
    mutate(if_connect = ifelse(n > n_spot_minimum | percent > 0.2, 1, 0)) %>%
    filter(if_connect == 1) 

# make graph
connected_region_graph = graph_from_data_frame(connected_df[,c('region1','region2')], directed = F)
# Get components (connected regions)
# https://stackoverflow.com/questions/56332850/how-can-i-extract-clusters-from-an-igraph-network
connected_region_graph_components = components(connected_region_graph)
threeD_cluster = connected_region_graph_components$membership

# Assign component to each region
connected_df = connected_df %>% mutate(
    volumn_cluster1 = connected_region_graph_components$membership[region1],
    volumn_cluster2 = connected_region_graph_components$membership[region2]
    )
write_tsv(connected_df, str_glue('{out_path}/4a_3D_cluster.tsv'))

# 4b. Assign to spatial plot
slice_library_name_match = sheet_multi$LibraryName %>% setNames(sheet_multi$PASTEFileName)
# Make spatial plot to see
spatial_map_df = 
    bind_rows(
        matching_df[,c('x_1','y_1','Filtered_tumor_regions_1','section1')] %>% rename(x=x_1, y=y_1, Filtered_tumor_regions=Filtered_tumor_regions_1, section=section1),
        matching_df[,c('x_2','y_2','Filtered_tumor_regions_2','section2')] %>% rename(x=x_2, y=y_2, Filtered_tumor_regions=Filtered_tumor_regions_2, section=section2)
        ) %>% 
        filter(
            Filtered_tumor_regions != '0'
            ) %>%
        mutate(section = factor(section, levels = mixedsort(unique(section)))) %>% 
        mutate(region = str_c(section,"-",Filtered_tumor_regions)) %>% 
        # Assign component
        mutate(volumn_cluster = threeD_cluster[match(region, names(threeD_cluster))]) %>%
        mutate(volumn_cluster = str_c('3D_',volumn_cluster)) %>% 
        # Add piece name
        mutate(PASTEpieceName = piece_use) %>%
        # add LibraryName
        mutate(LibraryName = slice_library_name_match[section])

write_tsv(spatial_map_df, str_glue('{out_path}/4b_3D_cluster_spatial_map.tsv'))

spatial_label_df = spatial_map_df %>% group_by(section,Filtered_tumor_regions) %>% summarize(across(where(is.numeric), mean))
p_spatial_3dcluster = ggplot(spatial_map_df, aes(x=x,y=-y, color=as.character(volumn_cluster))) +
    facet_wrap(~section) +
    geom_point(size = 0.5) + 
    geom_text(data = spatial_label_df, aes(label = Filtered_tumor_regions), color = 'gray20') + 
    cowplot::theme_cowplot() +
    coord_fixed() +
    theme(
        legend.position = "bottom"
        ) + 
    labs(title = '3D cluster')
ggsave(str_glue("{out_path}/4b_3D_cluster_spatial_map.pdf"), p_spatial_3dcluster, width = 10, height = 5, units = "in", dpi = 300)

# 5. save a simply 3D cluster label
volumn_cluster_clean_df = spatial_map_df %>%
    mutate(
        PASTEpieceName = piece_use
    ) %>% 
    dplyr::rename(
        PASTE_x = x,
        PASTE_y = y,
        PASTEFileName = section
    ) %>% select(LibraryName, Filtered_tumor_regions, volumn_cluster, PASTEFileName) %>% 
    distinct()

write_tsv(volumn_cluster_clean_df, str_glue('{out_path}/5_3D_cluster_label.tsv'))

#############################################################################################
#############################################################################################
# NEW. 2023/10/19 

# 6. Make sankey color by volume cluster --------------------------------------------------------
# Library
library(networkD3)
library(dplyr)
library(gtools)
# Assign color for each region based on volumn cluster
# fix NA. 
volumn_cluster_region_df = spatial_map_df %>% distinct(volumn_cluster, region) %>% 
    mutate(volumn_cluster = ifelse(is.na(volumn_cluster), 'NA', volumn_cluster)) 
# Get unique volumn cluster and assign color each 
color_volumn_cluster = unique(volumn_cluster_region_df$volumn_cluster) %>% 
    mixedsort %>% 
    setNames(colorspace::rainbow_hcl(length(.)),.)

# Assign color to each volumn cluster
my_color_volumn_group = define_colors(names(color_volumn_cluster), color_volumn_cluster)
group_color_df = volumn_cluster_region_df %>% setNames(c("node_name", "color_group"))
color_group = color_volumn_cluster
color_group_sankey_format = define_colors(names(color_group), color_group)

# Get links info
# NOTE. This format is important. AddDummyLinks uses these 5 columns
links_df = count_connection_df %>% mutate(value = n) %>% 
    ungroup %>% 
    select(section1, section2, region1, region2, value) %>%
    setNames(c("From_layer","To_layer",'From_node','To_node','value'))


### ------------------------------------------------------------
# A. Define layer order. Could be manually. If didn't provide or NULL, will use mixedsort to determine
layer_ranks = gtools::mixedsort(union(links_df$From_layer, links_df$To_layer))
links_df = AddDummyLinks(links_df) #, layer_ranks)
# ------------------------------------------------------------
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_df <- data.frame(
  node_name= union(as.character(links_df$From_node), as.character(links_df$To_node))
) 

# # Add color group to links, merged by region1
links_df = links_df %>% left_join(group_color_df, by = c('From_node' = 'node_name')) # link color determined by From_node
nodes_df = nodes_df %>% left_join(group_color_df, by = c('node_name' = 'node_name')) 
 
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links_df$IDsource <- match(links_df$From_node, nodes_df$node_name)-1 
links_df$IDtarget <- match(links_df$To_node, nodes_df$node_name)-1


# Make the Network
p <- sankeyNetwork(Links = links_df, Nodes = nodes_df,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "node_name",
              #iterations = 0 ,
              colourScale = color_group_sankey_format,
              LinkGroup = 'color_group',
              NodeGroup = 'color_group',
              sinksRight=FALSE)
p

# Save the result as html
require(htmlwidgets)
saveWidget(p, file=str_glue("{out_path}/6_sankey_color_group.html"))

require(webshot) # need conda install phantomjs
webshot(url = str_glue("{out_path}/6_sankey_color_group.html"), 
    file= str_glue("{out_path}/6_sankey_color_group.pdf"),
    vwidth = 700,
    vheight = 400,
)


# ----------------------------------------------------------------------------------------
# 7. Make sankey color by volume cluster --------------------------------------------------------
# with simpler node name
# Library
library(networkD3)
library(dplyr)
library(gtools)
# Assign color for each region based on volumn cluster
# fix NA. 
volumn_cluster_region_df = spatial_map_df %>% distinct(volumn_cluster, region) %>% 
    mutate(volumn_cluster = ifelse(is.na(volumn_cluster), 'NA', volumn_cluster)) 
# Get unique volumn cluster and assign color each 
color_volumn_cluster = unique(volumn_cluster_region_df$volumn_cluster) %>% 
    mixedsort %>% 
    setNames(colorspace::rainbow_hcl(length(.)),.)

# Assign color to each volumn cluster
my_color_volumn_group = define_colors(names(color_volumn_cluster), color_volumn_cluster)
group_color_df = volumn_cluster_region_df %>% setNames(c("node_name", "color_group"))
color_group = color_volumn_cluster
color_group_sankey_format = define_colors(names(color_group), color_group)

# Get links info
# NOTE. This format is important. AddDummyLinks uses these 5 columns
links_df = count_connection_df %>% mutate(value = n) %>% 
    ungroup %>% 
    select(section1, section2, region1, region2, value) %>%
    setNames(c("From_layer","To_layer",'From_node','To_node','value'))


### ------------------------------------------------------------
# A. Define layer order. Could be manually. If didn't provide or NULL, will use mixedsort to determine
layer_ranks = gtools::mixedsort(union(links_df$From_layer, links_df$To_layer))
links_df = AddDummyLinks(links_df) #, layer_ranks)
# ------------------------------------------------------------
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_df <- data.frame(
  node_name= union(as.character(links_df$From_node), as.character(links_df$To_node))
) 

# # Add color group to links, merged by region1
links_df = links_df %>% left_join(group_color_df, by = c('From_node' = 'node_name')) # link color determined by From_node
nodes_df = nodes_df %>% left_join(group_color_df, by = c('node_name' = 'node_name')) 
nodes_df = nodes_df %>% mutate(node_name_simple = str_extract(node_name, "-[0-9]+$") %>% str_remove('-')) # Add simply node name for just microregion
 
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links_df$IDsource <- match(links_df$From_node, nodes_df$node_name)-1 
links_df$IDtarget <- match(links_df$To_node, nodes_df$node_name)-1


# Make the Network
p <- sankeyNetwork(Links = links_df, Nodes = nodes_df,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "node_name_simple",
              #iterations = 0 ,
              colourScale = color_group_sankey_format,
              LinkGroup = 'color_group',
              NodeGroup = 'color_group',
              sinksRight=FALSE)
p

# Save the result as html
require(htmlwidgets)
saveWidget(p, file=str_glue("{out_path}/7_sankey_color_group_simple.html"))

require(webshot) # need conda install phantomjs
webshot(url = str_glue("{out_path}/7_sankey_color_group_simple.html"), 
    file= str_glue("{out_path}/7_sankey_color_group_simple.pdf"),
    vwidth = 700,
    vheight = 400,
)

