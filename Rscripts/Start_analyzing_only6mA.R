# 2025-05-19
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(dplyr)
library(IRanges)
library(ggplot2)
library(pheatmap)
library(rtracklayer)
# add some

# globals
project <- "Bio253"
dateToday <- format(Sys.Date(), "%Y%m%d")
desc <- "Prepare"

#
# data to read in
#


# read in the data
prX <- readRDS("../resources/SA6850_prXallWide_moreMetaInfo.rds") ### Please be aware that you may need to change the directory (i.e. ../resources/) if you have manually downloaded this repo
mtX <- readRDS("../resources/methylation_data_EGGnogAnnotated.rds")

# biocyc table here all genes are listed from left to rigth - no matter if + or - strand
biocycSA6850 <- read_tsv("../resources/Biocyc_allGenesSA6850.txt")
colnames(biocycSA6850) <- c("geneName", "Acc1", "start", "end", "product", "strand", "GOBP", "GOMF", "pw", "Acc2")

# what do we have in prX
colnames(prX)


biocycSA6850_slim <- biocycSA6850 %>%
    select(Acc2, start, end, strand) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))

# who is in the game what clones
table(mtX$group)

# Show densities of IPDratio and visualize different feature types
(pdfN <- paste(project,"_", desc, "IPDRatio_density_by_feature_", dateToday, ".pdf", sep = ""))
pdf(pdfN, width = 8, height = 6)
ggplot(mtX, aes(x = IPDRatio, fill = feature)) +
    geom_density(alpha = 0.5) +
    scale_x_log10() +
    labs(
        title = "Density of IPD Ratios by Feature Type",
        x = "IPD Ratio (log10 scale)",
        y = "Density",
        fill = "Feature Type"
    ) +
    theme_minimal()
dev.off()

# only 6mA are true valid, 5mC are FP calls
# filter for IPD ratio 2.8
# filter for 6mA



# globals
IPDrThreshold <- 2.8
featureOfInterest <- "m6A"

# join mtX 2 prX
# prX_wMtX <- left_join(x = prX, y = mtX, by = c("LocusTag" = "locus_tag"))
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))
table(mtX_wPrX$feature)

# filter the data for the specific clone and condition
mtX_wPrX_OfInterest <- mtX_wPrX %>%
    filter(feature == featureOfInterest) %>%
    filter(IPDRatio > IPDrThreshold)


table(mtX_wPrX_OfInterest$feature)
colnames(mtX_wPrX_OfInterest)
table(mtX_wPrX_OfInterest$group, mtX_wPrX_OfInterest$treatment)


# slim it down to better visualize
methylation_slim <- mtX_wPrX_OfInterest %>% select(start, strand, feature, motif, group, treatment) |> distinct()
table(methylation_slim$group, methylation_slim$treatment)


# Task0: show distribution along the chromosome

# visualize data as density on the chromosome with start as position
# Create the density plot with proper direction based on strand

# Histogram approach
ggplot() +
    geom_histogram(data = subset(methylation_slim, strand == "+"),
                   aes(x = start, y = after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    geom_histogram(data = subset(methylation_slim, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = paste0("Distribution of Methylation Sites by Strand (",")"),
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")



(fN <- paste(project,"_", desc,"_", "Distribution_mSites_alongChromosome", dateToday, ".pdf", sep = ""))
# save the plot
ggsave(fN, width = 20, height = 6)
# Conclusion ->  far from uniform distributed
# maybe investigate peaks and gaps further?


colnames(methylation_slim)
# Do binning of methylation sites to visualize as heatmap, we need all replicates
bin_size <- 100  # 10bp bins
methylation_binned_all <- mtX_wPrX_OfInterest %>% select(start, strand, feature, motif, group, treatment, replicate) %>%
    unite(condition, c(group, treatment, replicate), sep="_") %>%
    mutate(bin = floor(start / bin_size) * bin_size) %>%
    count(condition, strand, bin, name = "sites_in_bin") %>%
    complete(condition, strand, bin, fill = list(sites_in_bin = 0)) %>%
    pivot_wider(names_from = condition, values_from = sites_in_bin, values_fill = 0)


# make bin start, end
methylation_binned_all$start <- methylation_binned_all$bin
methylation_binned_all$endBin <- methylation_binned_all$bin + (bin_size -1) # end is start + bin_size - 1

# visualize the methylation bins with features
colnames(methylation_binned_all)
colnames(methylation_binned_all)[3:20]
mat <- methylation_binned_all[,3:20] # select only the methylation marks

rownames(mat) <- paste(methylation_binned_all$bin, methylation_binned_all$strand, sep = "_")
colnames(mat)

# Extract sample annotations
treatment <- gsub(".*_(.*)_.*", "\\1", colnames(mat))
clone <- gsub("_.*", "", colnames(mat))

annotation_col <- data.frame(
    Treatment = treatment,
    Clone = clone
)
rownames(annotation_col) <- colnames(mat)

# Define colors that exactly match your annotation values
ann_colors <- list(
    Treatment = c(
        TSB = "#d95f02",   # orange
        PASN = "#1b9e77"   # green
    ),
    Clone = c(
        `6850` = "#7570b3",   # purple
        SB0804 = "#e7298a",   # pink
        SB1002 = "#66a61e"    # green
    )
)

# Generate heatmap
pheatmap(
    mat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    main = "Methylation Binned Heatmap",
    fontsize_row = 6,
    fontsize_col = 10
)



# make a heatmap of the methylation marks
(pdfN <- paste(project,"_", desc, "_Heatmap_binned_", dateToday, ".pdf", sep = ""))
pdf(pdfN, 6,6)
pheatmap(mat,
         scale="none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = paste("Methylation Binned Heatmap"),
         fontsize_row = 6,
         fontsize_col = 10)

dev.off()

# calc and add column with difference to next position on the same strand
# mutate(diff_to_next_site = start - lag(start)) %>%

# identify peaks and gaps
# gaps!

methylation_data_with_diff <- methylation_slim %>%
    select(start, strand, feature, motif) %>%
    distinct() %>%
    arrange(start) %>%
    group_by(strand) %>%
    mutate(diff_to_next_site = lead(start) - start) %>%
    ungroup()

# inspect the diff values
# @NZ: do you think this is real?
# @NZ: What do you think could this mean?
table(table(methylation_data_with_diff$diff_to_next_site))

# visualize diff per strand as boxplot and use the same colors for strand
(pdfN <- paste(project,"_", desc, "_Boxplot_spacing_methylationSites_byStrand_", dateToday, ".pdf", sep = ""))
pdf(pdfN, width = 8, height = 6)
    ggplot(methylation_data_with_diff, aes(x = strand, y = diff_to_next_site, fill = strand)) +
    geom_boxplot() +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = paste0("Spacing between methylation sites"),
        x = "Strand",
        y = "Difference in Chromosome Position",
        fill = "Strand"
    ) +
    theme_minimal()
dev.off()

# peaks go for binned data!
# Create the density plot with proper direction based on strand
smoothin_param <- 0.04

ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = treatment),
                 alpha = 0.3, adjust = smoothin_param) +  # <-- less smoothing) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = treatment),
                 alpha = 0.3,
                 adjust = smoothin_param) +  # <-- less smoothing) +
    scale_fill_manual(values = c("PASN" = "darkblue", "TSB" = "yellow")) +
    facet_grid(group ~ .) +
    labs(
        title = "Density of Methylation Sites by Clone and Treatment",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")

(fN <- paste(project,"_", desc, "_Density_methylationSites_byClone_and_Treatment_", dateToday, ".pdf", sep = ""))
# save the plot
ggsave(fN, width = 20, height = 10)

# by treatment and clone
ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = group),
                 alpha = 0.3, adjust = smoothin_param) +  # <-- less smoothing) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = group),
                 alpha = 0.3,
                 adjust = smoothin_param) +  # <-- less smoothing) +
    scale_fill_manual(values = c("6850" = "lightblue", "SB0804" = "yellow", "SB1002" = "orange")) +
    facet_grid(treatment ~ .) +
    labs(
        title = "Density of Methylation Sites by Clone and Treatment",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")

(fN <- paste(project,"_", desc, "_Density_methylationSites_byTreatmentNClone_", dateToday, ".pdf", sep = ""))
# save the plot
ggsave(fN, width = 30, height = 10)


# do again binning without clones
# Do binning of methylation sites to visualize as heatmap, we need all replicates
bin_size <- 1000  # 100bp bins
methylation_binned_condition <- mtX_wPrX_OfInterest %>% select(start, strand, feature, motif, treatment) %>%
    mutate(bin = floor(start / bin_size) * bin_size) %>%
    count(treatment, strand, bin, name = "sites_in_bin") %>%
    complete(treatment, strand, bin, fill = list(sites_in_bin = 0)) %>%
    pivot_wider(names_from = treatment, values_from = sites_in_bin, values_fill = 0)


# find the peak region
methylation_binned_condition
methylation_binned_condition$TSBnPASN <- methylation_binned_condition$TSB + methylation_binned_condition$PASN
methylation_binned_condition$endBin <- methylation_binned_condition$bin + bin_size
methylation_binned_condition[order(methylation_binned_condition$TSBnPASN, decreasing = TRUE)[1:10],]



# loot into the gff to find features in bins?
mygff <- rtracklayer::readGFF("../resources/SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff", version = 3)
colnames(mygff)
mygffdf <- as.data.frame(mygff)
head(mygffdf)

# gffs are quite redundant -> slim it down
# we wanna keep one for proteins and one for genes

mygffdf_slim_gene <- mygffdf %>%
    filter(type %in% c("gene")) %>%
    select(type, start, end, strand, ID, Name, gene, locus_tag, protein_id) %>%
    distinct()

mygffdf_slim_protein <- mygffdf %>%
    filter(type %in% c("CDS")) %>%
    select(type, start, end, strand, ID, Name, gene, locus_tag, protein_id) %>%
    distinct()


# find genes in regions of peaks and gaps for StringDB enrichment?

# In peak regions

source("../Rscripts/HelperFunctions.R")
my_top_10_bins <- methylation_binned_condition[order(methylation_binned_condition$TSBnPASN, decreasing = TRUE)[1:10],]
features_in_top10_bins_protein <- find_overlapping_features(my_top_10_bins, mygffdf_slim_protein)
features_in_top10_bins_genes <- find_overlapping_features(my_top_10_bins, mygffdf_slim_gene)




# spacing regions of interest
SpacingCutoff <- 10000  # 10kb

methylation_withDiffs <- methylation_data_with_diff %>%
    arrange(desc(diff)) %>%
    mutate(hypoMethRegion_end = start + diff)

largest_diff_AllStrand_ofInterest <- methylation_withDiffs %>%
    filter(diff > SpacingCutoff) %>%
    select(start, hypoMethRegion_end, diff, strand, feature)


# Find genes overlapping with hypo-methylated regions
genes_in_hypo_regions <- find_overlapping_genes(biocycSA6850_slim, largest_diff_AllStrand_ofInterest)
genes_in_hypo_regions

# Some genes in biocyc do not have an RSAU_Number -> but NA --> follow up on these more?? what are these?
genes_in_hypo_regions <- genes_in_hypo_regions[!is.na(genes_in_hypo_regions)]


# Print results
cat("Number of Acc2 genes overlapping with hypo-methylated regions:", length(genes_in_hypo_regions), "\n")
print(genes_in_hypo_regions)

# Find the full details of these genes
overlapping_genes_details_all <- biocycSA6850_slim %>%
    filter(Acc2 %in% genes_in_hypo_regions)

# Print the details
print(overlapping_genes_details_all)

# bring in AGU
# here we parsed the description lines of the fasta file into a data frame
desc <- readRDS("../resources/descFasta.rds")
overlapping_genes_details_all <- left_join(overlapping_genes_details_all, desc, by = c("Acc2" = "LocusTag"))


# also join in proteomics expression data
# in this sheet we assembled the protein expression data w/ the different statistical analysis for different contrasts
prXdata <- read_tsv("../resources/SA6850_prXallWide_moreMeta.tsv")
genes_w_proteinExpression_in_hypoRegions <- left_join(x = overlapping_genes_details_all, y = prXdata, by = c("Acc2" = "LocusTag"))



# iBAQ value is a measurement for proteins in proteomics that can be used to
# compare different proteins to each other
# idea: compared to overall expression are the proteins often more or less abundant?
hist(log10(prXdata$mean_iBAQ))
hist(log10(genes_w_proteinExpression_in_hypoRegions$mean_iBAQ))

# combine these 2 histograms all iBAQs vs hypo iBAQs
ggplot() +
    geom_histogram(data = prXdata,
                   aes(x = log10(mean_iBAQ), y = after_stat(count), fill = "All Proteins"),
                   alpha = 0.5, bins = 20) +
    geom_histogram(data = genes_w_proteinExpression_in_hypoRegions,
                   aes(x = log10(mean_iBAQ), y = after_stat(count), fill = "Hypo Methylated Proteins"),
                   alpha = 0.5, bins = 20) +
    scale_fill_manual(values = c("All Proteins" = "blue", "Hypo Methylated Proteins" = "red")) +
    labs(
        title = paste0("Distribution of iBAQ values (log10) "),
        x = "log10(iBAQ)",
        y = "Count",
        fill = "Proteins"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")






# Look for Motif where 6mA is methylated for all OR for the ones in particular peaks?



# screen "back-ground" (DNA) for known (degenerated) motifs of methylation






