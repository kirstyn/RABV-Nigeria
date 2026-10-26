# ============================================================
# Pairwise genetic distance analysis
# ============================================================

library(ape)
library(ggplot2)
library(pheatmap)
library(patchwork)


# ============================================================
# Read alignment & sample metadata
# ============================================================

aln <- read.dna(
  "outputs/study_wgs.aln.fasta",
  format = "fasta"
)

meta <- read.csv("data/sample_accession_ids.csv")

# ============================================================
# Tidy names
# ============================================================

new_names <- gsub("/.*", "", labels(aln))
new_names[6] <- "AU105"
new_names <- meta$NCBI_accession[match(new_names, meta$Sample_id)]
rownames(aln) <- new_names

# ============================================================
# Calculate sequence coverage
# ============================================================

# Convert DNAbin alignment to a character matrix
aln_mat <- as.character(aln)

# Ensure it is a matrix with the same dimensions as the alignment
dim(aln_mat) <- dim(aln)

rownames(aln_mat) <- rownames(aln)
colnames(aln_mat) <- colnames(aln)

# Valid nucleotide positions
valid <- aln_mat %in% c(
  "a", "c", "g", "t",
  "A", "C", "G", "T"
)

# Make sure valid is a matrix
dim(valid) <- dim(aln_mat)

# Calculate percentage of genome with valid nucleotide data
coverage <- rowSums(valid) / ncol(aln_mat) * 100

# Add sequence names
names(coverage) <- rownames(aln_mat)

# Summary of sequence coverage
summary(coverage)

# View individual sequence coverage
coverage

# ============================================================
# Number of comparable sites for each pair
# ============================================================

# Pairwise number of positions where BOTH sequences
# have valid nucleotide data

n_sites <- matrix(
  NA,
  nrow = nrow(aln_mat),
  ncol = nrow(aln_mat),
  dimnames = list(
    rownames(aln_mat),
    rownames(aln_mat)
  )
)

for (i in seq_len(nrow(aln_mat))) {
  
  for (j in seq_len(nrow(aln_mat))) {
    
    n_sites[i, j] <- sum(
      valid[i, ] & valid[j, ]
    )
    
  }
}


# Summary of comparable sites across pairwise comparisons
summary(
  n_sites[upper.tri(n_sites)]
)


# Convert comparable sites to percentage of genome
n_sites_pct <- (
  n_sites / ncol(aln_mat)
) * 100

summary(
  n_sites_pct[upper.tri(n_sites_pct)]
)


# ============================================================
# Raw pairwise p-distance
# ============================================================

# pairwise.deletion = TRUE means that each pairwise
# comparison uses only sites where BOTH sequences have
# valid nucleotide data.

d <- dist.dna(
  aln,
  model = "raw",
  pairwise.deletion = TRUE
)


# Extract pairwise distances
dist_values <- as.vector(d)


# ============================================================
# Mean pairwise distance
# ============================================================

mean_dist <- mean(dist_values)

mean_dist
range(dist_values)


# ============================================================
# Calculate % identities
# ============================================================

# Mean pairwise nucleotide identity
mean_identity <- (1 - mean_dist) * 100

mean_identity


# Range/distribution of pairwise nucleotide identity
identity_values <- (1 - dist_values) * 100

summary(identity_values)


# Sensitivity analysis:
# sequences with >=80% genome coverage
# ============================================================

coverage_threshold <- 80

# Number of sequences retained
sum(coverage >= coverage_threshold)

# Alignment containing only high-coverage sequences
aln_80 <- aln[
  coverage >= coverage_threshold,
]


# Pairwise p-distance for high-coverage sequences
d_80 <- dist.dna(
  aln_80,
  model = "raw",
  pairwise.deletion = TRUE
)


# Extract distances
dist_values_80 <- as.vector(d_80)


# Mean distance for high-coverage sequences
mean_80 <- mean(dist_values_80)


# Mean identity for high-coverage sequences
mean_identity_80 <- (1 - mean_80) * 100


# ============================================================
# Compare all sequences vs high-coverage sequences
# ============================================================

comparison <- data.frame(
  Dataset = c(
    "All sequences",
    ">=80% genome coverage"
  ),
  N_sequences = c(
    nrow(aln),
    nrow(aln_80)
  ),
  Mean_p_distance = c(
    mean_dist,
    mean_80
  ),
  Mean_identity_percent = c(
    mean_identity,
    mean_identity_80
  )
)

comparison


# ============================================================
# Panel A: Histogram
# ============================================================

p1 <- ggplot(
  data.frame(distance = dist_values),
  aes(x = distance)
) +
  
  geom_histogram(
    bins = 30,
    colour = "black",
    fill = "grey80"
  ) +
  
  geom_vline(
    xintercept = mean_dist,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = mean_dist,
    y = Inf,
    label = paste0(
      "Mean = ",
      round(mean_dist, 4)
    ),
    vjust = 1.5,
    hjust = -0.1,
    size = 4
  ) +
  
  labs(
    x = "Pairwise p-distance",
    y = "Number of comparisons"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(),
    panel.grid = element_blank()
  )

p1


# ============================================================
# Panel B: Heatmap
# ============================================================

ph <- pheatmap(
  as.matrix(d),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  silent = TRUE
)

p2 <- wrap_elements(ph$gtable)

p2


# ============================================================
# Combine panels
# ============================================================

combined_plot <- p1 + p2 +
  
  plot_layout(
    widths = c(1, 1.5)
  ) +
  
  plot_annotation(
    tag_levels = "A"
  )

combined_plot


# ============================================================
# Export
# ============================================================

ggsave(
  "figures/pairwise_distance_figure.pdf",
  combined_plot,
  width = 12,
  height = 6
)

ggsave(
  "figures/pairwise_distance_figure.png",
  combined_plot,
  width = 12,
  height = 6,
  dpi = 600
)


library(adegenet)
library(here)

## Cluster analysis of genetic data, and lineage assignments

pat_matrix <- as.matrix(read.table(here("output/phylogenetic_distances/romblonSeq_24_patristicDist_matrix.txt"), check.names = FALSE))

# Plot frequency distribution of Hamming distances
jpeg(here("output/figures/Hamming_distance_frequencies.jpeg"))
hist(pat_matrix, 
     breaks = seq(0, max(pat_matrix) + 0.0005, by=0.0005),
     col = "cyan",
     xlab = "Hamming distance", main = "Distribution of frequencies")
dev.off()

# GENETIC CLUSTERS- ADEGENET
# The resulting cut-off of 0.0004 was selected by comparing resulting clusters with
# phylogenetic trees .
clust_all <- gengraph(d, show.graph=F, col.pal=rainbow, truenames=T, cex=0.3, cutoff=0.0004)
write.csv(clust_all$clust$membership, here("output/phylogenetic_distances/0.0004clust.csv"))

