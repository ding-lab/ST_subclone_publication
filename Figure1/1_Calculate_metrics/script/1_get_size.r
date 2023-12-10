# load sample and objects 
source('/PATH/TO/0_Load_All_Sample/script/load_samples.r')

# Parameters
analysis_date = ''
PATH_TO_ANALYSIS = ""
out_path = str_glue('{PATH_TO_ANALYSIS}/out/')
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# Get size of each region
# scaling_factor_for_between_spot_space = 2.1 # This is proven wrong. Stop using this! 
spot_size = (55/2)^2 * pi
per_spot_area_v116 = 50^2*sqrt(3)/2*1/2 * 8 # 8 small triangle with 50um edge. + 2 extra tranigle same size
# -- convertion ratio v1, per spot area: 1 / per_spot_area (um^2 / spot) * 10^6 = 10^6 = 1 / 8660 * 10^6 = 115.4
sample_size_df = imap(st_list, function(st_use, sample_id){
    st_use@meta.data %>% count(genetic_clone, Filtered_tumor_regions) %>%
        mutate(sample_name = sample_id) %>% 
        dplyr::rename(n_spot = n) %>% 
        mutate(
            region_size_mm2_v116 = n_spot *  per_spot_area_v116 / 10^6
            )
}) %>% bind_rows()

# add annotation
sample_size_df = sample_size_df %>% mutate(TumorTME = ifelse(Filtered_tumor_regions == '0', 'TME', 'Tumor'))

# Add size group
# Remake size group based on new criteria - V20231024
high_cutoff = 250
low_cutoff = 25
sample_size_df = sample_size_df %>% mutate(size_group = case_when(
  n_spot <= low_cutoff ~ 'small',
  n_spot <= high_cutoff ~ 'mid',
  TRUE ~ 'large'
))
sample_size_df = sample_size_df %>% mutate(high_cutoff = high_cutoff, low_cutoff= low_cutoff)

# write out
write_tsv(sample_size_df, file.path(out_path, str_glue('ST_microregion_size_v{analysis_date}.tsv')))


########################################
# Count tumor percent, density, overall and for each tumor size group
# 1. Tumor percent =  n_Tumor / n_Tumor + n_TME
# 2. Overall microregion density = n_microregion / total area (in spots) * convertion_ratio (in spot /mm^2)
# -- convertion ratio v1, per spot area: 1 / per_spot_area (um^2 / spot) * 10^6 = 1 / 8660 * 10^6 = 115.4
# 3. microregion density for each tumor size group:
# -- small: n_microregion / total area (in spots) * convertion_ratio (in spot /mm^2)
# -- mid: n_microregion / total area (in spots) * convertion_ratio (in spot /mm^2)
# -- large: n_microregion / total area (in spots) * convertion_ratio (in spot /mm^2)

# v116
# Convertion ratio 115.5 (n_spot / mm^2) or 8660 (um^2 / spot)
convertion_ratio_use =  115.5
ST_microregion_density_convertion_116 = sample_size_df %>% 
    group_by(sample_name) %>% 
    summarize(
        tumor_percent = sum(n_spot[Filtered_tumor_regions != '0']) / sum(n_spot),
        n_micoregion_overall = n_distinct(Filtered_tumor_regions[Filtered_tumor_regions != '0']),
        n_micoregion_small = n_distinct(Filtered_tumor_regions[Filtered_tumor_regions != '0' & size_group == 'small']),
        n_micoregion_mid = n_distinct(Filtered_tumor_regions[Filtered_tumor_regions != '0' & size_group == 'mid']),
        n_micoregion_large = n_distinct(Filtered_tumor_regions[Filtered_tumor_regions != '0' & size_group == 'large']),
        overall_density = n_micoregion_overall / sum(n_spot) * convertion_ratio_use, # in mm2
        small_density = n_micoregion_small / sum(n_spot) * convertion_ratio_use,
        mid_density = n_micoregion_mid / sum(n_spot) * convertion_ratio_use,
        large_density = n_micoregion_large / sum(n_spot) * convertion_ratio_use,
    ) %>% 
    # Add overall density per slide
    mutate(across(ends_with('density'), ~`*`(.,6.5^2), .names = "{.col}_per_slide"))
    
 write_tsv(ST_microregion_density_convertion_116, file.path(out_path, '2_ST_microregion_density_convertion_ratio116.tsv'))

########################################
# PLOT and Check result
library(patchwork)
# plot Estimated n micoregion vs observed n micoregion
## v116
density_df_use = ST_microregion_density_convertion_116
# label those difference greater than 20%
density_df_use = density_df_use %>% mutate(label_sample = ifelse(abs(overall_density_per_slide - n_micoregion_overall) / n_micoregion_overall > 0.5, sample_name, NA))
p_convert116 = map(c('overall','small','mid','large'), function(plot_group){
    df_use = density_df_use %>% select(sample_name, tumor_percent,label_sample, contains(plot_group))
    ggplot(df_use, aes(
        x = .data[[str_glue('n_micoregion_{plot_group}')]], 
        y = .data[[str_glue('{plot_group}_density_per_slide')]]
        )) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = label_sample), size = 2) +
    geom_smooth(method = 'lm') +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    labs(x = 'Observed n microregion', y = 'Estimated n microregion', 
        title = plot_group
        ) +
    theme_bw()
}) %>% wrap_plots()

pdf(file.path(out_path, '2_ST_microregion_density_convertion_ratio116.pdf'), width = 8, height = 8)
p_convert116 + plot_annotation(
        title = 'Estimated vs Observed n microregion using convertion ratio 116.3',
        subtitle = "Sample with estimated n > 50% of observed n are labeled"
        )
dev.off()

#######################################
# 20231030 - 
# From people feedback ... use per mm^2 version
## uses 116 spot/mm^2 version
ST_microregion_density_convertion_116_annotated = 
    left_join(ST_microregion_density_convertion_116, 
        sheet_use %>% select(LibraryName, sample_type_tumor, cancer_type),
              by = c('sample_name' = 'LibraryName')) 

sample_order = ST_microregion_density_convertion_116_annotated %>% arrange(overall_density_per_slide) %>% pull(sample_name)

### Also make one for per slide
ST_microregion_density_convertion_116_long_per_mm = 
    ST_microregion_density_convertion_116_annotated %>% 
    select(!contains('overall')) %>%
    pivot_longer(cols = ends_with('_density'), names_to = 'density_type', values_to = 'density')
    
p_density_116_perslide_bar = ST_microregion_density_convertion_116_long_per_mm %>%
    mutate(sample_type_tumor = factor(sample_type_tumor, levels = c('Primary', 'Metastasis'))) %>%
ggplot(.,
    aes(y = factor(sample_name, levels = sample_order), x = density * 100, fill = density_type)) +
    geom_bar(stat = 'identity') + 
    geom_line() +
    #scale_y_continuous(labels = scales::comma) +
    facet_grid(cancer_type + sample_type_tumor~., space= 'free', scales = 'free') +
    theme_bw() + 
    scale_fill_manual(values = col_size_group_plus) + 
    # remove minor grid
    theme(panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = 'bottom') +
    labs(x = 'Density(n microregion / cm2)', y = '',
        title = 'Estimated n microregion density per slide using convertion ratio 112.87')

pdf(file.path(out_path, '3d_ST_microregion_density_convertion_ratio116_long_per_cm2_stackedbar.pdf'), width = 5, height = 6)
p_density_116_perslide_bar
dev.off()
