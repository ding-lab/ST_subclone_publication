
library(patchwork)

# External functions
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script_git/Clustering/function_reoderbyhcluster.r')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/ST_subclone/38-GenomicSubclone/5_GSEA/script/src/function_TestAndFixAddModuleScore.r')
## ---------- Plot Function ---------- ##
# 1. Dotplot/Heatmap with annotation
AnnoDotPlot = function(obj, group.by = 'seurat_clusters', features, annotation_idents, label_ident = T, 
    mode = c("Dot","Heatmap"),
    cluster_row = T,
    cluster_col = T,
    # title
    title = NULL, subtitle = NULL,
    # Module score
    AddModuleScore = T, ModuleScoreHeight = 2,
    # Highlight based on value cutoff
    highlight_tiles = F,
    highlight_cutoff_quantile = 0.75, 
    highlight_color = '#333333',
    highlight_thickness = 1,
    # Column to splot
    split_column = NULL,
    ...){
    message("annotation_idents:", annotation_idents, "Take values from @meta.data")
    # A. Main Dot/Heatmap plot
    pdata = DotPlot(obj, group.by = group.by, features = features) %>% .$data
    # A0. Filter out no expression gene
    pdata = pdata %>% filter(!is.nan(avg.exp.scaled))
    pdata = pdata %>% filter(!is.na(features.plot)) # This is weird need check
    # A1. Hierarchical clustering
    if(cluster_row) pdata = pdata %>% ReorderByHCluster(ident_column = 'id', groupby_column = 'features.plot', value_column = 'avg.exp.scaled')
    if(cluster_col) pdata = pdata %>% ReorderByHCluster(ident_column = 'features.plot', groupby_column = 'id', value_column = 'avg.exp.scaled')
    # test = pdata %>% ReorderSplitByHCluster(ident_column = 'id', groupby_column = 'features.plot', value_column = 'avg.exp.scaled', split_by_vector = )
    # return(test) 
    
    p_dot_exp = pdata %>% ggplot(aes(x = id, y = features.plot)) 
    # A2. Dot of Heatmap plot
    mode = match.arg(mode)
    if(mode == 'Dot'){
        p_dot_exp = p_dot_exp + 
            geom_point(aes(color = avg.exp.scaled, size = pct.exp)) + 
            scale_color_gradient2(low = '#3333DD', mid = '#E0E0E0', high = '#DD3333', midpoint = 0) 
    }else if(mode == 'Heatmap'){
        #Heatmap
        p_dot_exp = p_dot_exp +
        geom_tile(aes(fill = avg.exp.scaled)) +
        scale_fill_gradient2(low = '#3333DD', mid = '#E0E0E0', high = '#DD3333', midpoint = 0) 
    }
    # Highlight specific tiles
    if(highlight_tiles){
        highlight_cutoff_values = quantile(pdata$avg.exp.scaled, probs = highlight_cutoff_quantile)
        pdata_highlight = pdata %>% filter(avg.exp.scaled > highlight_cutoff_values)
        p_dot_exp = p_dot_exp +
            geom_tile(data = pdata_highlight, aes(fill = avg.exp.scaled), color = highlight_color, linejoin= "round", linewidth = highlight_thickness) 
    }
    # Theme Adjustments
    p_dot_exp = p_dot_exp + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
    # Get groupby orders
    order_group_by = levels(pdata$id)
    # Get final feature count
    p_dot_exp_nrow = length(unique(pdata$features.plot))


    # B. Create column annotation bars
    # Use idetns in meta.data
    annotation_df = FetchData(obj, vars = c(group.by, annotation_idents))
    annotation_collapsed_df = annotation_df %>% 
        group_by(.data[[group.by]]) %>% 
        summarize(across(all_of(annotation_idents), ~paste(sort(unique(.)), collapse = ' '))) %>% 
        mutate({{group.by}} := factor(.data[[group.by]], levels = order_group_by)) # Reoder groupby
    # # Create bar/column plot
    p_bar_nrow = length(unique(annotation_idents))
    p_bar_list = map(annotation_idents, function(anno_ident){
        p_bar = annotation_collapsed_df[, c(group.by, anno_ident)] %>% 
            ggplot(aes(x = .data[[group.by]], y = anno_ident, fill = .data[[anno_ident]])) + 
            geom_tile() +
            theme_void() + 
            #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  
            # Remove x axis 
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
            # Add y text back
            theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5)) + 
            # Use color scale colorspace::rainbow_hcl(n)
            scale_fill_manual(values = colorspace::rainbow_hcl(n = length(unique(annotation_collapsed_df[[anno_ident]]))))
        if(label_ident) p_bar = p_bar + geom_text(aes(label = str_wrap(.data[[anno_ident]], width = 4))) 
        return(p_bar)
    }) %>% setNames(annotation_idents)

    # B1. Make current height arragment
    plot_height_arrangement = c(rep(1,p_bar_nrow), p_dot_exp_nrow)

    # B. Add module score
    if(AddModuleScore){
        # First test if need to fix object
        obj = TestAndFixSeuratForAddModuleScore(obj)
        # Plot
        p_modulescore = ModuleScoreBoxplot(obj, group.by = group.by, features_plt = features) + 
            theme_bw() + 
            # Remove x axis
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
            # Remove minor grid
            theme(panel.grid.minor = element_blank()) +
            # Add y text back
            theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
        # Reorder column
        p_modulescore$data = p_modulescore$data %>% 
            mutate({{group.by}} := factor(.data[[group.by]], levels = order_group_by)) # Reoder groupby
        # update plot arrangement
        p_bar_list = c(list(ModuleScore = p_modulescore), p_bar_list)
        plot_height_arrangement = c(ModuleScoreHeight, plot_height_arrangement) # Append height of module score to top 
    }

    # C. Combined
    p_all = wrap_plots(c(p_bar_list, list(Dotplot= p_dot_exp)), ncol = 1, heights = plot_height_arrangement, guides = "collect") &
         # put annotation to bottom
        theme(legend.position = 'bottom')  
    # D, Titles
    p_all = p_all + plot_annotation(title = title, subtitle = subtitle, 
        theme = theme(
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            plot.subtitle = element_text(hjust = 0.5, face = 'italic'))
        )
    return(p_all)
}

# Module score boxplot
ModuleScoreBoxplot = function(obj, group.by = 'seurat_clusters', features_plt){
    #browser()
    message("Calculating Module score")
    obj_tmp = AddModuleScore(obj, features = list(ModuleScore=features_plt))
    obj_tmp@meta.data[['ModuleScore']] = obj_tmp@meta.data[['Cluster1']]
    
    # Module score boxplot
    FetchData(obj_tmp, vars = c(group.by, 'ModuleScore')) %>% 
        mutate({{group.by}} := as.character(.data[[group.by]])) %>%
        ggplot(aes(x = .data[[group.by]], y = ModuleScore, fill = .data[[group.by]], group = .data[[group.by]])) +
        geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
