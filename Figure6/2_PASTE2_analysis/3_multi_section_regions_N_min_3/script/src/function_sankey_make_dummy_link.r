# This helper function is used to make dummy links for sankey plot so that nodes are at the right level
# The format for layer node rank df is:
#    Layer    Node        Layer_rank
#    <chr>    <chr>            <dbl>
#  1 Th1H3U1  Th1H3U1-4            1
#  2 Th1H3U1  Th1H3U1-3            1
#  3 Th1H3U1  Th1H3U1-6            1
#  4 Th1H3U12 Th1H3U12-6           3
#  5 Th1H3U12 Th1H3U12-10          3
#  6 Th1H3U12 Th1H3U12-5           3
#  7 Th1H3U2  Th1H3U2-2            2
#  8 Th1H3U2  Th1H3U2-1            2
#  9 Th1H3U2  Th1H3U2-5            2

make_dummy_link_one_node = function(missing_node, layer_node_rank_df){
    # 0. get layer ranks
    layer_ranked = layer_node_rank_df %>% distinct(Layer,Layer_rank) %>% arrange(Layer_rank) %>% pull(Layer)
    # 1. get the rank
    missing_node_rank = layer_node_rank_df %>% filter(Node == missing_node) %>% pull(Layer_rank) 
    missing_node_layer = layer_node_rank_df %>% filter(Node == missing_node) %>% pull(Layer)
    # First, create n-1 to actual node with missing layers. e.g rank = 4, dummy_layer_3 to missing_node
    dummy_links_df = data.frame(
            From_layer = layer_ranked[[missing_node_rank-1]], 
            To_layer = missing_node_layer, 
            From_node = str_glue("Dummy_node_{missing_node_rank-1}"),
            To_node = missing_node,
            value = 1
        )
    # 2nd, if rank > 2, Create all n-1 prior dummy layer links. e.g. rank = 4, create 1-2, 2-3, 2 layers. 
    if(missing_node_rank > 2){
        dummy_links_prior_df = map(seq(1, missing_node_rank-2, 1), function(missing_rank){
        data.frame(
            From_layer = layer_ranked[[missing_rank]], 
            To_layer = layer_ranked[[missing_rank+1]], 
            From_node = str_glue("Dummy_node_{missing_rank}"),
            To_node = str_glue("Dummy_node_{missing_rank+1}"),
            value = 1
            )
        }) %>% bind_rows() %>% distinct()
        dummy_links_df = bind_rows(dummy_links_prior_df, dummy_links_df)
    }
    
    return(dummy_links_df)
}


# Main function to add dummy links
AddDummyLinks = function(links, layer_ranks=NULL){
    # A. set parameter
    if(is.null(layer_ranks)){
        message("layer_ranks not provided. automatically determined by mixedsort")
        layer_ranks = gtools::mixedsort(union(links$From_layer, links$To_layer))
        message("layer_ranks: ", toString(layer_ranks))
    }
    first_layer_name = layer_ranks[[1]]
    # function to help sort: https://stackoverflow.com/questions/32378108/using-gtoolsmixedsort-or-alternatives-with-dplyrarrange
    mixedrank = function(x) order(gtools::mixedorder(x))

    # B. create Layer-Node-LayerRank map to help create dummy layer
    layer_node_rank_df = links %>% {bind_rows(
        select(., contains('From_')) %>% setNames(c("Layer", "Node")),
        select(., contains('To_'))  %>% setNames(c("Layer", "Node"))
        )} %>% distinct %>% arrange(mixedrank(Layer), Node) %>% 
        mutate(Layer_rank = factor(Layer, levels = layer_ranks) %>% as.numeric) 


    # C. Find nodes missing levels 
    true_first_layer_nodes = links %>% filter(From_layer == first_layer_name) %>% pull(From_node)
    current_top_level_nodes = setdiff(links$From_node, links$To_node ) 
    nodes_missing_levels = setdiff(current_top_level_nodes, true_first_layer_nodes)

    # D. create all dummry rows for all missing nodes
    dummy_links_df = map(nodes_missing_levels, function(missing_node){
        make_dummy_link_one_node(missing_node, layer_node_rank_df)
    }) %>% bind_rows() %>% distinct() 
    message("Adding follow dummy links:")
    print(dummy_links_df)
    # Add to links
    links = bind_rows(links, dummy_links_df) %>% distinct()
    return(links)
}
