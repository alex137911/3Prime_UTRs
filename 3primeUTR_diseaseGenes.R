# Code to compare 3' UTR composition across disease gene groups

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(showtext))

# Fonts
font_add_google("Poppins", family = "Poppins")
showtext_auto()
showtext_opts(dpi = 100)

# -------------------------------------------------------------------
# IMPORT DATA
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/3C/3Prime_UTRs/Data/preliminaryAnalysis.R"
inDir  <- sprintf("%s/preliminaryAnalysis.R", dirname(inpath))

# Import processed 3' UTR exon length data (from preliminaryAnalysis.R)
threeUTR_exons <- read_delim(file = sprintf("%s/threeUTR_exonLengths.tsv", inDir), delim = "\t")
threeUTR_exons <- as.data.frame(threeUTR_exons)
message(sprintf("%i transcripts", nrow(threeUTR_exons)))

inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/3C/3Prime_UTRs/Data/Input"
inDir  <- sprintf("%s/Input", dirname(inpath))

# Developmental disorder dominant/recessive genes
# Downloaded on August 17, 2023 (https://www.ebi.ac.uk/gene2phenotype/downloads)
ddGenes <- read_csv(file.path(inpath, "DDG2P_17_8_2023.csv.gz"))

# COSMIC Cancer Gene Census
# Downloaded on August 17, 2023 (https://cancer.sanger.ac.uk/cosmic/download)
cancerGenes <- read_tsv(file.path(inpath, "Cosmic_CancerGeneCensus_v98_GRCh38.tsv.gz"))

# Dosage sensitive data
# R. L. Collins et al., “A cross-disorder dosage sensitivity map of the human genome,”
# Cell, vol. 185, no. 16, pp. 3041-3055.e25, Aug. 2022, doi: 10.1016/j.cell.2022.06.036.
# Downloaded August 17, 2023 (https://zenodo.org/record/6347673)
doseSensitive <- read_tsv(file.path(inpath, "Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"))

# Read in NHS PanelApp green genes
greenGenes_NHS <- read_table(file = sprintf("%s/greenGenes_NHS.txt", inDir), col_names = FALSE)
greenGenes_NHS <- greenGenes_NHS %>% gather(value = "gene_name") %>% select(-key)

# -------------------------------------------------------------------
# Output directory
outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/3C/3Prime_UTRs/Data"
outDir <- sprintf("%s/Data/3primeUTR_diseaseGenes", dirname(outpath))
if(!file.exists(outDir)) dir.create(outDir)

write.table(ddGenes, file = sprintf("%s/developmental_disorderGenes.tsv", outDir),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -------------------------------------------------------------------
# Merge 3' UTR exon with NHS green genes
threeUTR_exonMerged <- threeUTR_exons %>%
  mutate(green_gene = gene_name %in% greenGenes_NHS$gene_name)

# ---------------------------------
# Merge cancer genes
cancerGenes <- cancerGenes %>% rename(gene_name = GENE_SYMBOL)
  
# Mutation abbreviations: https://cancer.sanger.ac.uk/cosmic/help/census#abbrev
cancerGenes_filtered <- dplyr::filter(cancerGenes, grepl("N|F|Mis|D", MUTATION_TYPES))

# Filter for oncogenes (includes genes that are classified as oncogenes + TSG)
oncogene <- dplyr::filter(cancerGenes_filtered, grepl("oncogene", ROLE_IN_CANCER))

# Filter for tumour suppressor genes (without oncogenes)
TSG <- dplyr::filter(cancerGenes_filtered, grepl("TSG", ROLE_IN_CANCER))
TSG_filtered <- TSG %>% anti_join(oncogene, by = "gene_name")

# Merge oncogenes
threeUTR_exonMerged <- threeUTR_exonMerged %>%
  mutate(oncogene = gene_name %in% oncogene$gene_name)

# Merge TSGs
threeUTR_exonMerged <- threeUTR_exonMerged %>%
  mutate(TSG = gene_name %in% TSG_filtered$gene_name)

# ---------------------------------
# Merge developmental disorder genes
ddGenes <- ddGenes %>% rename(gene_name = "gene symbol")
ddGenes <- ddGenes %>% rename(allelic_requirement = "allelic requirement")
ddGenes <- ddGenes[ !ddGenes$`confidence category` == "limited", ] # 368 genes

# Dominant DD genes
dominantDD <- dplyr::filter(ddGenes, grepl("monoallelic", allelic_requirement))

# Recessive DD genes
recessiveDD <- dplyr::filter(ddGenes, grepl("biallelic", allelic_requirement))

# Merge dominant
threeUTR_exonMerged <- threeUTR_exonMerged %>%
  mutate(dominant_dd = gene_name %in% dominantDD$gene_name)

# Merge recessive
threeUTR_exonMerged <- threeUTR_exonMerged %>%
  mutate(recessive_dd = gene_name %in% recessiveDD$gene_name)

# ---------------------------------
# Merge dosage sensitive genes
doseSensitive <- doseSensitive %>% rename(gene_name = "#gene")

# Haploinsufficient genes
haploinsufficient <- doseSensitive[ doseSensitive$pHaplo >= 0.86, ] # 2987 genes

# Triplosensitive genes
triplosensitive   <- doseSensitive[ doseSensitive$pTriplo >= 0.94, ] # 1559 genes

# Merge haploinsufficient
threeUTR_exonMerged <- threeUTR_exonMerged %>%
  mutate(haploinsufficient = gene_name %in% haploinsufficient$gene_name)

# Merge triplosensitive
threeUTR_exonMerged <- threeUTR_exonMerged %>%
  mutate(triplosensitive = gene_name %in% triplosensitive$gene_name)

# -------------------------------------------------------------------
# Boxplot
boxPlot_colours <- c("#BAE4C4", "#B7E3D4", "#A4DCD7", "#89D3D3", "#89C1D3", 
                     "#4EB3D3", "#2B8CBE", "#0868AC", "#094081", "#042042")

# Separate data frames for each disease gene group
df_all <- threeUTR_exonMerged

df_green_genes <- filter(threeUTR_exonMerged, green_gene == TRUE)

df_cancer <- threeUTR_exonMerged %>% 
  filter(oncogene == TRUE | TSG == TRUE) %>%
  mutate(category = ifelse(oncogene == TRUE, "Oncogene", "TSG"))

df_dev_disorder <- threeUTR_exonMerged %>% 
  filter(recessive_dd == TRUE | dominant_dd == TRUE) %>%
  mutate(category = ifelse(recessive_dd == TRUE, "Recessive DD", "Dominant DD"))

df_dosage_sensitive <- threeUTR_exonMerged %>% 
  filter(haploinsufficient == TRUE | triplosensitive == TRUE) %>%
  mutate(category = ifelse(haploinsufficient == TRUE, "Haploinsufficient", "Triplosensitive"))

# Combine all data frames
df_combined <- bind_rows(
  df_all %>% mutate(plot_category = "All Genes"),
  df_green_genes %>% mutate(plot_category = "NHS Green Genes"),
  df_cancer %>% mutate(plot_category = "Cancer"),
  df_dev_disorder %>% mutate(plot_category = "Developmental Disorder"),
  df_dosage_sensitive %>% mutate(plot_category = "Dosage Sensitive")
)

# Plot
# Calculate the median of df_all
median_val <- median(df_all$total_length)

# Wrap label 
df_combined$plot_category[df_combined$plot_category == "NHS Green Genes"] <- "NHS Green\nGenes"

# Define custom y-axis labels
table(df_combined$plot_category)

df_combined <- df_combined %>%
  mutate(category = case_when(
    plot_category == "All Genes" ~ "All Genes (n=18872)",
    plot_category == "Cancer" & category == "TSG" ~ "TSG (n=191)",
    plot_category == "Cancer" & category == "Oncogene" ~ "Oncogene (n=163)",
    plot_category == "Developmental Disorder" & category == "Recessive DD" ~ "Recessive (n=1242)",
    plot_category == "Developmental Disorder" & category == "Dominant DD" ~ "Dominant (n=851)",
    plot_category == "Dosage Sensitive" & category == "Haploinsufficient" ~ "Haploinsufficient (n=2858)",
    plot_category == "Dosage Sensitive" & category == "Triplosensitive" ~ "Triplosensitive (n=1490)",
    plot_category == "NHS Green\nGenes" & is.na(category) ~ "Green Genes (n=3508)",
    TRUE ~ category
))

ggplot(df_combined, aes(y = category, x = total_length, fill = category)) +
  geom_boxplot() +
  geom_vline(aes(xintercept = median_val), linetype = "dotted", color = "black", size = 0.8) +
  facet_grid(plot_category ~ ., scales = "free_y", space = "free_y",) +
  xlab("3' UTR Length (bp)") + ylab("Disease Gene Groups") +
  coord_cartesian(xlim = c(0, 5000), expand = TRUE) +
  scale_fill_manual(values = boxPlot_colours) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(family = "Poppins"),
        axis.text = element_text(family = "Poppins"),
        plot.title = element_text(family = "Poppins"),
        plot.margin = margin(50, 50, 20, 20),  
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)))

showtext_opts(dpi = 600)

ggsave("3' UTR Length Across Disease Genes.png", plot = last_plot(), device = "png", 
       width = 19.0, height = 8.42, dpi = 600)