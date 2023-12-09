# Cancer color from Clara:
source('/diskmnt/Projects/Users/cliu/pancan_ST/color.R')

color_cancer_discover_full = c(Breast = "#d386b7", Colorectal = "#b75420", Pancreas = "#5d9db9"
)
color_cancer_discover = c(BRCA = "#d386b7", CRC = "#b75420", PDAC = "#5d9db9")
color_cancer_all = c(color_cancer_discover,
    ccRCC = "#BB9946",
    RCC = "#BB9946", # "#22bc6c",
    CHOL = "#bfa8db", #Cholangiocarcinoma = "#bfa8db",
    UCEC = "#f39c12"
)

color_n_section = c('5' = '#3862a5','4'='#588fcc','3'='#85b1dd','2'='#abcfee','1'='#dfe9f1')

color_organ = c(
    "Pancreas" = "#8EA0C7",
    "Kidney" = "#E9BE55", 
    "Breast" = "#d386b7",
    "Liver" = "#55ae90",
    "Lung" = "#88bbee", 
    "Colon" = "#E9916E", 
    "Bile Duct" = "#b15e22",
    "Uretus" = "#f39c12",
    "Cholangio" = "#bfa8db",
    "Prostate" = "#2980b9"
)

color_tumor_type = c(
    Metastasis = "#C54F69",
    Primary = "#55749D"
)

color_assay = c('FFPE' = '#865C8F', 'OCT' = '#99AE62')

color_cohort = c('Validation' = "#0072B2", 'Discovery' = '#D55E00')

# A distinct color reference set 
# by GPT and manual selection
color_distinct = c(
"#e74c3c", #(bright red)
"#c0392b", #(deep red)
"#b15e22", #(burnt orange)
"#b9471f", #(rusty orange)
"#e7b3a3", #(pale peach)
"#e89c84", #(coral)
"#f39c12", #(orange)
"#f1c40f", #(bright yellow)
"#e5b15f", #(mustard yellow)
"#e4bf80", #(pale gold)
"#d3b98d", #(pale yellow-brown)
"#d9c9aa", #(pale beige)
"#27ae60", #(deep green)
"#1abc9c",
"#a8bfb1", #(gray-green)
"#8299a2", #(slate blue)
"#c0d2c5", #(pale sage green)
"#5fa4ad", #(seafoam green)
"#a5b1c2", #(pale blue-gray)
"#70b6c1", #(pale teal)
"#2980b9", #(bright blue)
"#bfa8db", #(pale lavender)
"#b894d0", #(pale lilac)
"#d8a8d9", #(pale pink)
"#d1a6bc", #(pale mauve)
"#c88fbf", #(mauve-pink)
"#d68ab5" #(soft pink)
)

color_size_group = c('small' = '#4E6F99', 'mid' = '#65D196', 'large' = '#E7E431')
color_primary_met = c('Primary' = '#7B9AD0', 'Metastasis' = '#EE9933')


# c('Brain' = '#0072B2', 'Kidney' = '#D55E00', 'Liver' = '#CC79A7', 'Lung' = "#CCff12",'Lymph Node' = '#11ffbb',"Muscle" = "#22FF13",'Skin' = '#5572B2', 'Spleen' = '#D55E33', 'Stomach' = '#CC7912'),