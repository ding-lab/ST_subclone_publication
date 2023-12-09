# Add average expression of genetic subclone 
CalculateCloneExpressionLevel = function(obj, features, assay = 'SCT', group.by = 'genetic_clone'){
    message("Note, currently use 0 as TME and clone to detect tumor columns")
    message("Assay = ", assay)
    # Calculate average expression of genetic subclone
    exp_df = obj %>% 
        AverageExpression(
            assays = assay, slot = 'data', 
            group.by =  group.by,
            features = features) %>% 
            .[[assay]] %>% 
        as.data.frame %>%
        rownames_to_column('Gene')
    # Add Min, max and difference 
    exp_df %>% 
        rowwise() %>% 
        # Get Min and Max
        mutate(
            max_tumor_value = max(c_across(contains('clone'))),
            max_tumor = unlist(pmap(across(contains('clone')), ~names(c(...)[which.max(c(...))]))),
            min_tumor_value = min(c_across(contains('clone'))),
            min_tumor = unlist(pmap(across(contains('clone')), ~names(c(...)[which.min(c(...))])))
            # ^^ https://stackoverflow.com/questions/17735859/for-each-row-return-the-column-name-of-the-largest-value
            ) %>% 
        # Add difference
        mutate(
            min_max_tumor_ratio = max_tumor_value / min_tumor_value,
            max_tme_ratio = max_tumor_value / `0`
        )
}