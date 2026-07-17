#------------------------------
# Libraries
#------------------------------

library(ape)
library(treeio)
library(ggtree)
library(tidyverse)
library(ggstar)


#------------------------------
# 1. Import tree
#------------------------------

tree <- read.tree(
  "data/fig3_sequences/combined_set_sequences.wOutgp.aln.fasta.contree"
)
tree <- drop.tip(
  tree,
  tip = "DQ837463"
)

#------------------------------
# 2. Root tree by outgroup
#------------------------------

outgroup <- "EU293119.1"

tree_rooted <- root(
  tree,
  outgroup = outgroup,
  resolve.root = TRUE
)

tree_no_outgroup <- drop.tip(
  tree_rooted,
  "EU293119.1"
)


#------------------------------
# 3. Import metadata
#------------------------------

metadata <- read_delim(
  "data/fig3_sequences/combined_set_metadata.txt",
  delim = "\t"
)


#------------------------------
# 4. Link metadata to tree
#------------------------------

p <- ggtree(tree_no_outgroup)

p$data <- p$data %>%
  left_join(
    metadata,
    by = c("label" = "locus")
  )


# check unmatched tips

p$data %>%
  filter(isTip, is.na(country)) %>%
  select(label)


#------------------------------
# 5. Plot tree
#------------------------------
p$data <- p$data %>%
  mutate(
    tip_labels = str_c(
      label,
      country,
      sep = "_"
    )
  )
p$data <- p$data %>%
  mutate(
    study = str_detect(authors, regex("Ukande", ignore_case = TRUE)) &
      str_detect(authors, regex("Kia", ignore_case = TRUE))
  )

clade_nodes <- p$data %>%
  filter(isTip) %>%
  group_by(EPA_major_clade) %>%
  summarise(
    node = if (n() > 1) {
      getMRCA(tree_no_outgroup, label)
    } else {
      first(node)
    },
    n = n(),
    .groups = "drop"
  )
clade_nodes <- clade_nodes %>%
  mutate(
    clade_label = recode(
      EPA_major_clade,
      "Bats" = "Bats DR"
    )
  )

clade_nodes_nocosmo <- clade_nodes%>%
  filter(EPA_major_clade != "Cosmopolitan")


minor_clade_nodes <- p$data %>%
  filter(isTip) %>%
  group_by(EPA_minor_clade) %>%
  summarise(
    node = if (n() > 1) {
      getMRCA(tree_no_outgroup, label)
    } else {
      first(node)
    },
    n = n(),
    .groups = "drop")%>% 
  filter(
    EPA_minor_clade != "DR",
    !is.na(EPA_minor_clade)
  )

p$data <- p$data %>%
  
  mutate(
    
    support = as.numeric(na_if(label, ""))
    
  )

p_tree <- p +
  geom_tree(
    linewidth = 0.5
  ) +
  geom_point2(
    aes(
      subset = !isTip & support > 80
    ),
    shape = 23,        
    size = 2,
    alpha=0.6,
    fill = "grey60",
    colour = "black"
  )+
  geom_tippoint(
    aes(fill = EPA_major_clade),
    shape = 21,
    size = 2,
    colour = "black"
  ) +
  geom_tiplab(
    aes(label = tip_labels),
    size = 4,
    align = FALSE, 
    offset=0.005
  ) +
  geom_star(
    data = function(x) subset(x, study == TRUE),
    aes(y = y),
    starshape = 1,
    size = 4,
    fill = "gold",
    colour = "black"
  )+
  geom_treescale(
    x = 0.002,
    y = 10,
    width = 0.02
  ) +
  geom_cladelab(node = clade_nodes_nocosmo$node, 
                label = clade_nodes_nocosmo$clade_label,  
                offset = 0.09,
                fontface = "bold")+
  geom_cladelab(node = minor_clade_nodes$node, 
                label = minor_clade_nodes$EPA_minor_clade,  
                offset = 0.09,
                fontface = "bold.italic")+
  # geom_text2(
  #   aes(
  #     subset = node %in% clade_nodes$node,
  #     label = support
  #   ),
  #   data = left_join(
  #     p$data,
  #     clade_nodes %>% select(node, clade_label),
  #     by = "node"
  #   ),
  #   size = 3,
  #   fontface = "bold.italic"
  # )+
  theme_tree() +
  theme(
    legend.position = c(0.15, 0.8),   
    legend.background = element_rect(
      fill = "white",
      colour = NA
    ),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.25, "cm")
  );p_tree

p_tree <- p_tree +
  theme(
    plot.margin = margin(
      t = 10,
      r = 50,
      b = 10,
      l = 20
    )
  ) +
  coord_cartesian(clip = "off")
p_tree <- p_tree +
  labs(
    fill = "Major clade"
  )
# PDF (vector; preferred for journals)
ggsave(
  "figures/RABV_tree_clades.pdf",
  plot = p_tree,
  width = 12,
  height = 8,
  units = "in",
  device = cairo_pdf
)

# TIFF (high-resolution raster)
ggsave(
  "figures/RABV_tree_clades.tiff",
  plot = p_tree,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# PNG (optional)
ggsave(
  "figures/RABV_tree_clades.png",
  plot = p_tree,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

