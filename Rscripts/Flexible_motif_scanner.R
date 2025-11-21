#!/usr/bin/env Rscript
# Flexible Motif Scanner with Degenerate Base Support
# Searches genome for motifs with IUPAC ambiguity codes

# load libraries
library(Biostrings)
library(dplyr)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# functions

# Function to find all matches of a degenerate motif
find_degenerate_motif <- function(genome_seq, pattern) {

    # Convert pattern to DNAString with IUPAC codes
    pattern_dna <- DNAString(pattern)

    # Search on + strand
    cat("  Searching + strand...\n")
    plus_matches <- matchPattern(pattern_dna, genome_seq,
                                 with.indels = FALSE,
                                 fixed = FALSE)  # Allow IUPAC ambiguity
    plus_positions <- start(plus_matches)
    plus_sequences <- as.character(plus_matches)

    # Search on - strand (reverse complement of pattern)
    cat("  Searching - strand...\n")
    pattern_rc <- reverseComplement(pattern_dna)
    minus_matches <- matchPattern(pattern_rc, genome_seq,
                                  with.indels = FALSE,
                                  fixed = FALSE)
    minus_positions <- start(minus_matches)
    minus_sequences <- as.character(minus_matches)

    # Create results data frame
    results <- data.frame(
        position = c(plus_positions, minus_positions),
        strand = c(rep("+", length(plus_positions)),
                   rep("-", length(minus_positions))),
        sequence = c(plus_sequences, minus_sequences),
        stringsAsFactors = FALSE
    )

    # Sort by position
    results <- results[order(results$position), ]

    return(list(
        results = results,
        pattern_rc = as.character(pattern_rc)
    ))
}


#
# start of the analysis
#
#


# globals
project <- "Bio253"
dateToday <- format(Sys.Date(), "%Y%m%d")
desc <- "getBackgroundDensity"



# get mtx
MtX <- readRDS("../resources/methylation_data.rds")

# globals
IPDrThreshold <- 2.8
featureOfInterest <- "m6A"

# filter the data for the specific clone and condition
mtX_OI <- MtX %>%
    filter(feature == featureOfInterest) %>%
    filter(IPDRatio > IPDrThreshold)

# slim it down to better visualize
methylation_slim <- mtX_OI %>% select(start, strand) |> distinct()


# Load methylation positions
meth_data <- methylation_slim

# Load genome sequence
genome <- readDNAStringSet("../resources/SA_6850_ncb i_sequence.fasta")
genome_seq <- genome[[1]]  # Extract the DNAString object
genome_length <- length(genome_seq)
genome_length





# AI: claude.ai
# ============================================================================
# 1. LOAD DATA
# ============================================================================



# CHANGE THIS TO YOUR DESIRED MOTIF:
motif_pattern <- "RVANNNNNNNTBY"
# motif_pattern <- "ANNNNNNNT"
# motif_pattern <- "RNANNNNNNTNC"

# For your second analysis, change it to:
# motif_pattern <- "RVANNNNNNTHC"
# Which translates to: [G|A][G|A|C]ANNNNNNNT[T|G|C][C|T]

# Position within motif where methylation occurs (1-based)
# For ANNNNNNT, if methylation is on the A, set to 1
# If it's after 4 bases, set to 4, etc.
# methylation_offset <- 3  # Adjust based on your motif
methylation_offset <- 3  # Adjust based on your motif

# Matching tolerance for position matching (bp)
position_tolerance <- 0

cat("===================================================\n")
cat("FLEXIBLE MOTIF GENOME SCANNER\n")
cat("===================================================\n\n")
cat("Motif pattern:", motif_pattern, "\n")
cat("Motif length:", nchar(motif_pattern), "bp\n")
cat("Methylation position:", methylation_offset, "\n\n")

# ============================================================================
# 2. LOAD DATA
# ============================================================================

cat("Genome length:", genome_length, "bp\n")
cat("Detected methylation sites:", nrow(meth_data), "\n")
cat("  + strand:", sum(meth_data$strand == "+"), "\n")
cat("  - strand:", sum(meth_data$strand == "-"), "\n\n")

# ============================================================================
# 3. SEARCH FOR MOTIF WITH IUPAC CODES
# ============================================================================

cat("Scanning genome for motif pattern...\n")

# Function to find all matches of a degenerate motif
find_degenerate_motif <- function(genome_seq, pattern) {

    # Convert pattern to DNAString with IUPAC codes
    pattern_dna <- DNAString(pattern)

    # Search on + strand
    cat("  Searching + strand...\n")
    plus_matches <- matchPattern(pattern_dna, genome_seq,
                                 with.indels = FALSE,
                                 fixed = FALSE)  # Allow IUPAC ambiguity
    plus_positions <- start(plus_matches)
    plus_sequences <- as.character(plus_matches)

    # Search on - strand (reverse complement of pattern)
    cat("  Searching - strand...\n")
    pattern_rc <- reverseComplement(pattern_dna)
    minus_matches <- matchPattern(pattern_rc, genome_seq,
                                  with.indels = FALSE,
                                  fixed = FALSE)
    minus_positions <- start(minus_matches)
    minus_sequences <- as.character(minus_matches)

    # Create results data frame
    results <- data.frame(
        position = c(plus_positions, minus_positions),
        strand = c(rep("+", length(plus_positions)),
                   rep("-", length(minus_positions))),
        sequence = c(plus_sequences, minus_sequences),
        stringsAsFactors = FALSE
    )

    # Sort by position
    results <- results[order(results$position), ]

    return(list(
        results = results,
        pattern_rc = as.character(pattern_rc)
    ))
}

# Perform the search
search_results <- find_degenerate_motif(genome_seq, motif_pattern)
all_motifs <- search_results$results
pattern_rc <- search_results$pattern_rc

cat("\n=== SEARCH RESULTS ===\n")
cat("Total motif occurrences:", nrow(all_motifs), "\n")
cat("  + strand:", sum(all_motifs$strand == "+"), "\n")
cat("  - strand:", sum(all_motifs$strand == "-"), "\n")
cat("\nReverse complement pattern:", pattern_rc, "\n\n")

# Show first few matches as examples
if (nrow(all_motifs) > 0) {
    cat("First 10 matches:\n")
    print(head(all_motifs, 10))
    cat("\n")
}

# ============================================================================
# 4. CHECK METHYLATION STATUS
# ============================================================================

cat("Checking methylation status for each motif occurrence...\n\n")

# Calculate expected methylation position for each motif
# For + strand: position + offset - 1
# For - strand: position is still the start of the motif on the genome,
# but we need to account for the motif being read in reverse
all_motifs$expected_meth_pos <- ifelse(
    all_motifs$strand == "+",
    all_motifs$position + methylation_offset - 1,
    # For minus strand, the methylation offset is counted from the END of the motif
    # because the motif is read in reverse on that strand
    all_motifs$position + nchar(motif_pattern) - methylation_offset
)

# Initialize methylation status
all_motifs$is_methylated <- FALSE
all_motifs$closest_meth_dist <- NA

# Check each motif occurrence
for (i in 1:nrow(all_motifs)) {
    expected_pos <- all_motifs$expected_meth_pos[i]
    strand <- all_motifs$strand[i]

    # Find methylation sites on same strand
    strand_meth <- meth_data[meth_data$strand == strand, ]

    if (nrow(strand_meth) > 0) {
        # Calculate distances to all methylation sites
        distances <- abs(strand_meth$start - expected_pos)
        min_dist <- min(distances)

        all_motifs$closest_meth_dist[i] <- min_dist

        # Mark as methylated if within tolerance
        if (min_dist <= position_tolerance) {
            all_motifs$is_methylated[i] <- TRUE
        }
    }
}

# ============================================================================
# 5. CALCULATE STATISTICS
# ============================================================================

cat("=== METHYLATION COVERAGE STATISTICS ===\n\n")

total_motifs <- nrow(all_motifs)
methylated_motifs <- sum(all_motifs$is_methylated)
unmethylated_motifs <- total_motifs - methylated_motifs

if (total_motifs > 0) {
    coverage <- (methylated_motifs / total_motifs) * 100
} else {
    coverage <- 0
    cat("WARNING: No motifs found! Check your pattern.\n\n")
}

cat("Motif pattern:", motif_pattern, "\n")
cat("Total occurrences:", total_motifs, "\n")
cat("Methylated:", methylated_motifs, "\n")
cat("Unmethylated:", unmethylated_motifs, "\n")
cat("Coverage:", round(coverage, 2), "%\n\n")

# Strand-specific statistics
if (total_motifs > 0) {
    stats_by_strand <- all_motifs %>%
        group_by(strand) %>%
        summarise(
            total = n(),
            methylated = sum(is_methylated),
            unmethylated = sum(!is_methylated),
            coverage_pct = round((methylated / total) * 100, 2)
        )

    cat("Strand-specific statistics:\n")
    print(as.data.frame(stats_by_strand))
    cat("\n")
}

# Sequence variant analysis
if (total_motifs > 0) {
    cat("=== SEQUENCE VARIANTS ===\n")
    cat("Most common sequences matching the pattern:\n")
    seq_counts <- sort(table(all_motifs$sequence), decreasing = TRUE)
    print(head(seq_counts, 20))
    cat("\n")

    # Check methylation rate for each variant
    cat("Methylation rate by sequence variant (top 10):\n")
    variant_meth <- all_motifs %>%
        group_by(sequence) %>%
        summarise(
            count = n(),
            methylated = sum(is_methylated),
            meth_rate = round((methylated / count) * 100, 1)
        ) %>%
        arrange(desc(count))

    print(head(as.data.frame(variant_meth), 10))
    cat("\n")
}

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

if (total_motifs > 0 && methylated_motifs > 0) {

    # Plot 1: Methylation coverage by strand
    p1 <- ggplot(all_motifs, aes(x = strand, fill = is_methylated)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#27AE60"),
                          labels = c("Unmethylated", "Methylated"),
                          name = "Status") +
        labs(title = paste0("Methylation Coverage: ", motif_pattern),
             x = "Strand",
             y = "Proportion") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
        scale_y_continuous(labels = scales::percent) +
        geom_text(stat = 'count', aes(group = is_methylated,
                                      label = after_stat(count)),
                  position = position_fill(vjust = 0.5))

    fN <- paste(motif_pattern, "_", project,"_",desc, "_methylation_coverage_by_strand.png", sep="")
    ggsave(fN, p1, width = 8, height = 6, dpi = 300)
    print(p1)

    # Plot 2: Genome-wide distribution of motifs (your favorite!)
    p2 <- ggplot(all_motifs, aes(x = position, y = 0, color = is_methylated)) +
        geom_point(alpha = 0.6, size = 2) +
        facet_wrap(~strand, ncol = 1) +
        scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#27AE60"),
                           labels = c("Unmethylated", "Methylated"),
                           name = "Status") +
        labs(title = paste0("Genome-wide Distribution: ", motif_pattern),
             x = "Genomic Position (bp)",
             y = "") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) +
        scale_x_continuous(labels = scales::comma)

    fN <- paste(motif_pattern, "_", project,"_",desc, "_motif_distribution_genome.png", sep="")
    ggsave(fN, p2, width = 14, height = 6, dpi = 300)
    print(p2)

    # Plot 3: Count comparison
    count_data <- data.frame(
        Category = c("Total Motifs", "Methylated", "Unmethylated"),
        Count = c(total_motifs, methylated_motifs, unmethylated_motifs),
        Color = c("Total", "Methylated", "Unmethylated")
    )

    p3 <- ggplot(count_data, aes(x = Category, y = Count, fill = Color)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Count), vjust = -0.5, size = 6, fontface = "bold") +
        scale_fill_manual(values = c("Total" = "#3498DB",
                                     "Methylated" = "#27AE60",
                                     "Unmethylated" = "#E74C3C")) +
        labs(title = paste0("Motif Methylation Status: ", motif_pattern),
             subtitle = paste0("Coverage: ", round(coverage, 1), "%"),
             x = "",
             y = "Count") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              legend.position = "none") +
        ylim(0, max(count_data$Count) * 1.15)

    fN <- paste(motif_pattern, "_", project,"_",desc, "_methylation_counts.png", sep="")
    ggsave(fN, p3, width = 8, height = 6, dpi = 300)
    print(p3)

    # Plot 4: Top sequence variants with methylation status
    if (nrow(variant_meth) > 0) {
        top_variants <- head(variant_meth, 15)

        p4 <- ggplot(top_variants, aes(x = reorder(sequence, count), y = count)) +
            geom_bar(stat = "identity", aes(fill = meth_rate)) +
            geom_text(aes(label = paste0(meth_rate, "%")), hjust = -0.2, size = 3) +
            coord_flip() +
            scale_fill_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#27AE60",
                                 midpoint = 50, name = "Methylation\nRate (%)") +
            labs(title = paste0("Top Sequence Variants: ", motif_pattern),
                 x = "Sequence",
                 y = "Count") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

        fN <- paste(motif_pattern, "_", project,"_",desc, "_sequence_variants.png", sep="")
        ggsave(fN, p4, width = 10, height = 8, dpi = 300)
        print(p4)
    }
}

# ============================================================================
# 6B. DENSITY OVERLAY PLOTS
# ============================================================================

if (total_motifs > 0 && methylated_motifs > 0) {
    cat("\n=== GENERATING DENSITY OVERLAY PLOTS ===\n")

    # Prepare data for density plots
    # Methylation marks
    meth_plus <- meth_data %>% filter(strand == "+")
    meth_minus <- meth_data %>% filter(strand == "-")

    # Motif occurrences
    motif_plus <- all_motifs %>% filter(strand == "+")
    motif_minus <- all_motifs %>% filter(strand == "-")

    # Calculate bin width for histogram (adjust this for your genome size)
    # For a ~3Mb genome, 5kb bins give good resolution
    bin_width <- 2000

    # Create density overlay for + strand using histograms for high resolution
    # + Strand Density Overlay
    p5 <- ggplot() +
        geom_density(data = motif_plus,
                     aes(x = position, fill = "Motif occurrences"),
                     alpha = 0.4, adjust = 0.001) +
        geom_density(data = meth_plus,
                     aes(x = start, fill = "Methylation marks"),
                     alpha = 0.6, adjust = 0.001) +
        scale_fill_manual(values = c("Motif occurrences" = "#3498DB",
                                     "Methylation marks" = "#27AE60"),
                          name = "Feature") +
        labs(title = "Density Distribution: + Strand",
             subtitle = paste0("Motif: ", motif_pattern, " (density overlay, adjust=1.2)"),
             x = "Genomic Position (bp)", y = "Density") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              plot.subtitle = element_text(hjust = 0.5, size = 11),
              legend.position = "top") +
        scale_x_continuous(labels = scales::comma)

    fN <- paste(motif_pattern, "_", project,"_",desc, "_density_overlay_plus_strand.png", sep="")
    ggsave(fN, p5, width = 14, height = 6, dpi = 300)
    print(p5)

    # Create density overlay for - strand using histograms
    # - Strand Density Overlay
    p6 <- ggplot() +
        geom_density(data = motif_minus,
                     aes(x = position, fill = "Motif occurrences"),
                     alpha = 0.4, adjust = 0.001) +
        geom_density(data = meth_minus,
                     aes(x = start, fill = "Methylation marks"),
                     alpha = 0.6, adjust = 0.001) +
        scale_fill_manual(values = c("Motif occurrences" = "#3498DB",
                                     "Methylation marks" = "#E74C3C"),
                          name = "Feature") +
        labs(title = "Density Distribution: - Strand",
             subtitle = paste0("Motif: ", motif_pattern, " (density overlay, adjust=1.2)"),
             x = "Genomic Position (bp)", y = "Density") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              plot.subtitle = element_text(hjust = 0.5, size = 11),
              legend.position = "top") +
        scale_x_continuous(labels = scales::comma)

    fN <- paste(motif_pattern, "_", project,"_",desc, "_density_overlay_minus_strand.png", sep="")
    ggsave(fN, p6, width = 14, height = 6, dpi = 300)
    print(p6)

    # Combined Density Overlay Both Strands
    motif_all_labeled <- all_motifs %>%
        mutate(type = "Motif occurrences")
    meth_all_labeled <- meth_data %>%
        rename(position = start) %>%
        mutate(type = "Methylation marks")

    # Combined density plot (both strands) using histograms
    # Prepare combined data
    motif_all_labeled <- all_motifs %>%
        mutate(type = "Motif occurrences")
    meth_all_labeled <- meth_data %>%
        rename(position = start) %>%
        mutate(type = "Methylation marks")

    p7 <- ggplot() +
        geom_density(data = motif_all_labeled,
                     aes(x = position, fill = type),
                     alpha = 0.4, adjust = 0.001) +
        geom_density(data = meth_all_labeled,
                     aes(x = position, fill = type),
                     alpha = 0.6, adjust = 0.001) +
        facet_wrap(~strand, ncol = 1) +
        scale_fill_manual(values = c("Motif occurrences" = "#3498DB",
                                     "Methylation marks" = "#9B59B6"),
                          name = "Feature") +
        labs(title = "Density Distribution: Both Strands",
             subtitle = paste0("Motif: ", motif_pattern, " (density overlay, adjust=1.2)"),
             x = "Genomic Position (bp)", y = "Density") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              plot.subtitle = element_text(hjust = 0.5, size = 11),
              legend.position = "top") +
        scale_x_continuous(labels = scales::comma)

    fN <- paste(motif_pattern, "_", project,"_",desc, "_density_overlay_both_strands.png", sep="")
    ggsave(fN, p7, width = 14, height = 8, dpi = 300)
    print(p7)

    cat("Density overlay plots generated!\n")
}

# ============================================================================
# 7. EXPORT RESULTS
# ============================================================================

# Save all results
if (total_motifs > 0) {
    fN <- paste(motif_pattern, "_", project,"_",desc, "_motif_methylation_status.txt", sep="")
    write.table(all_motifs,
                fN,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)

    # Save unmethylated motifs
    if (unmethylated_motifs > 0) {
        unmethylated <- all_motifs[!all_motifs$is_methylated, ]
        fN <- paste(motif_pattern, "_", project,"_",desc, "_unmethylated_motifs.txt", sep="")
        write.table(unmethylated,
                    fN,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE)
    }

    # Save methylated motifs
    if (methylated_motifs > 0) {
        methylated <- all_motifs[all_motifs$is_methylated, ]
        fN <- paste(motif_pattern, "_", project,"_",desc, "_methylated_motifs.txt", sep="")
        write.table(methylated,
                    fN,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE)
    }

    # Save sequence variant analysis
    if (exists("variant_meth")) {
        fN <- paste(motif_pattern, "_", project,"_",desc, "_sequence_variant_analysis.txt", sep="")
        write.table(variant_meth,
                    fN,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE)
    }
}

# Save summary report
fN <- paste(motif_pattern, "_", project,"_",desc, "_motif_scan_summary.txt", sep="")
sink(fN)
cat("===================================================\n")
cat("FLEXIBLE MOTIF GENOME SCAN SUMMARY\n")
cat("===================================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Motif pattern:", motif_pattern, "\n")
cat("Pattern length:", nchar(motif_pattern), "bp\n")
cat("Reverse complement:", pattern_rc, "\n")
cat("Methylation offset:", methylation_offset, "\n")
cat("Position tolerance:", position_tolerance, "bp\n\n")
cat("Genome length:", format(genome_length, big.mark = ","), "bp\n\n")
cat("=== SEARCH RESULTS ===\n")
cat("Total occurrences:", total_motifs, "\n")
if (total_motifs > 0) {
    cat("Motif density:", round(total_motifs / (genome_length / 1000), 2), "per kb\n")
    cat("Average spacing:", round(genome_length / total_motifs), "bp\n\n")
}
cat("=== METHYLATION STATISTICS ===\n")
cat("Methylated motifs:", methylated_motifs, "\n")
cat("Unmethylated motifs:", unmethylated_motifs, "\n")
cat("Overall coverage:", round(coverage, 2), "%\n\n")
if (total_motifs > 0) {
    cat("=== STRAND STATISTICS ===\n")
    print(as.data.frame(stats_by_strand))
    cat("\n")
    cat("=== TOP SEQUENCE VARIANTS ===\n")
    print(head(as.data.frame(variant_meth), 15))
}
sink()

cat("\n===================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("===================================================\n")
if (total_motifs > 0) {
    cat("Output files generated:\n")
    cat("  - motif_methylation_status.txt\n")
    cat("  - methylated_motifs.txt\n")
    if (unmethylated_motifs > 0) {
        cat("  - unmethylated_motifs.txt\n")
    }
    cat("  - sequence_variants_analysis.txt\n")
    cat("  - methylation_coverage_by_strand.png\n")
    cat("  - motif_distribution_genome.png (your favorite!)\n")
    cat("  - methylation_counts.png\n")
    cat("  - sequence_variants.png\n")
    cat("  - density_overlay_plus_strand.png\n")
    cat("  - density_overlay_minus_strand.png\n")
    cat("  - density_overlay_both_strands.png\n")
    cat("  - motif_scan_summary.txt\n")
} else {
    cat("No motifs found matching pattern:", motif_pattern, "\n")
    cat("Try adjusting your pattern or check your input files.\n")
}
cat("===================================================\n")

