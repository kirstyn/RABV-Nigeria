library(ape)
library(ggplot2)
library(pheatmap)


aln <- read.dna(
  "outputs/study_wgs.aln.fasta",
  format = "fasta"
)

# raw p-distance

d <- dist.dna(
  aln,
  model = "raw",
  pairwise.deletion = TRUE
)

# mean pairwise distance

mean_dist <- mean(as.matrix(d)[upper.tri(as.matrix(d))])

mean_dist

# Extract upper triangle only

dist_values <- as.vector(d)

# Histogram
ggplot(data.frame(distance = dist_values),
       aes(x = distance)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = mean(dist_values),
             linetype = "dashed") +
  labs(
    x = "Pairwise p-distance",
    y = "Number of comparisons",
    subtitle = paste("Mean distance =", round(mean(dist_values), 4))
  ) +
  theme_bw()

library(pheatmap)

pheatmap(
  as.matrix(d),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Pairwise genetic distance matrix"
)
