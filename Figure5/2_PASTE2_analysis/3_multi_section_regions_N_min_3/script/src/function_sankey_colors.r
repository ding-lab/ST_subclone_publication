# usage: 
# my_color = define_colors(c('Th1H3U1','Th1H3U2'), c('#1f77b4','#ff7f0e'))
# Then in Sankey
# p <- sankeyNetwork(Links = links, Nodes = nodes,
#               Source = "IDsource", Target = "IDtarget",
#               Value = "value", NodeID = "name",
#               #iterations = 0 ,
#               colourScale = my_color,
#               sinksRight=FALSE)

define_colors = function(ident, colors){
    ident_string = ident %>%  paste(collapse = "','")
    colors_string = colors %>%  paste(collapse = "','")
    my_color <- str_glue("d3.scaleOrdinal() .domain(['{ident_string}']) .range(['{colors_string}'])")
    return(my_color)
}