
library(ape)
library(treeio)
library(ggtree)
library(tidyverse)
library(phangorn)
library(pals)
library(colorspace)
library(ggstar)
library(patchwork)
library(cowplot)
library(ggtreeExtra)


#------------------------------
# 1. Import tree
#------------------------------

tree <- read.tree("outputs/Africa2_allSeq_tree/RABV_AF2_iqtree.contree")

# midpoint root
tree <- midpoint(tree)


#------------------------------
# 2. Import metadata
#------------------------------

metadata <- read_csv("data/8jul26_all_afr2_metadata.csv") %>%
  mutate(
    country = str_replace(
      country,
      "^C\\?\\?te d'Ivoire$",
      "Cote d'Ivoire"
    )
  )



#------------------------------
# 4. Create ggtree object + metadata
#------------------------------

p <- ggtree(tree)

p$data <- p$data %>%
  left_join(
    metadata,
    by = c("label" = "locus")
  ) %>%
  left_join(
    node_data,
    by = "node"
  )


#------------------------------
# 5. Colours
#------------------------------

# country colours
n_country <- length(unique(p$data$country))

country_cols <- setNames(
  desaturate(
    glasbey(n_country),
    amount = 0.2
  ),
  unique(p$data$country)
)

colScale <- scale_fill_manual(
  name = "Country",
  values = country_cols,
  guide = guide_legend(
    ncol = 5,
    override.aes = list(
      size = 2.5
    )
  )
)

saveRDS(country_cols, "R/phylogeny-country-cols.rds")

#------------------------------
# 6. Identify study sequences
#------------------------------

p$data <- p$data %>%
  mutate(
    study = str_detect(authors, regex("Ukande", ignore_case = TRUE)) &
      str_detect(authors, regex("Kia", ignore_case = TRUE))
  )



#------------------------------
# 7. Plot- full view
#------------------------------

p_full <- p +
  geom_tree(
    linewidth = 0.4
  ) +
  geom_hilight(
    node = mrca_node,
    fill = "grey",
    alpha = 0.3

  ) +
  # 
  # geom_cladelabel(
  #   node = mrca_node,
  #   label = "Study sequences",
  #   offset = 0.01
  # )+
  geom_tippoint(
    aes(fill = country),
    shape = 21,
    size = 1.5,
    colour = "grey30",
    alpha = 0.8
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
    y = 800,
    width = 0.005
  ) +
  theme_tree() +
  colScale +
  theme(legend.position = "bottom")+
  labs(fill = "Country"); p_full

p_full <- p_full +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 3)
    )
  )

# plot- zoomed in 
study_tips <- p$data %>%
  filter(isTip, study == TRUE) %>%
  pull(label)

mrca_node <- getMRCA(tree, study_tips)

mrca_node

p_zoom <- viewClade(p, node = mrca_node)

p_zoom$data <- p_zoom$data %>%
  mutate(
    geo_loc_clean = geo_loc %>%
      str_remove(" State") %>%
      na_if("NA"),
    tip_label = if_else(
      is.na(geo_loc_clean),
      label,
      paste(label, geo_loc_clean, sep = ":")
    )
  )

p_zoom_plot <- p_zoom +
  geom_tree(linewidth = 0.5) +
  # geom_hilight(
  #   node = mrca_node,
  #   fill = "grey",
  #   alpha = 0.3
  #   
  # )+
  geom_point2(
    aes(subset = !isTip & as.numeric(label) > 80),
    size = 1.5,
    fill = "grey40",
    shape=23
  ) +
  geom_tippoint(
    aes(fill = country),
    shape = 21,
    size = 2.5,
    colour = "grey30"
  ) +
  geom_star(
    data = function(x) subset(x, study == TRUE),
    aes(y = y),
    starshape = 1,
    size = 4,
    fill = "gold",
    colour = "black"
  ) +
  geom_tiplab(
    aes(label = tip_label),
    size = 2.5,
    offset = 0.0001
  )+
  geom_treescale(
    x = 0.03,
    y = 891,
    width = 0.005
  )+
  theme_tree() +
  theme(legend.position = "none")+
  scale_fill_manual(values = country_cols);p_zoom_plot

p_zoom_plot2<- p_zoom_plot +
  coord_cartesian(
    xlim = c(0.026, 0.06),
    ylim = c(860, 892.5),
    clip = "off"
  );p_zoom_plot2

# fix margins a bit
p_full <- p_full +
  theme(
    plot.margin = margin(
      t = 30,   # top
      r = 20,
      b = 20,
      l = 20
    )
  ) +
  coord_cartesian(clip = "off")

p_full <- p_full +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.25, "cm"),
    legend.spacing.x = unit(0.05, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      nrow = 3,
      byrow = TRUE,
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  );p_full

p_full <- p_full +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title.position = "left",
    # reduce space around legend items
    legend.key.size = unit(0, "cm"),
    legend.key.width = unit(0, "cm"),
    legend.key.height = unit(0, "cm"),
    # reduce gaps between items
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    # reduce space between title and items
    legend.title.align = 0,
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_legend(
      nrow = 4,
      byrow = TRUE,
      keywidth = unit(0, "cm"),
      keyheight = unit(0, "cm"),
      label.position = "right",
      label.hjust = 0,
      override.aes = list(
        size = 3
      )
    )
  );p_full

    
p_zoom_plot_highlight <- p_zoom_plot2 + geom_hilight(
  node = mrca_node,
  fill = "grey",
  alpha = 0.3
)+
  guides(fill = "none") +
  theme(legend.position = "none")

p_full + p_zoom_plot_highlight +
  plot_layout(
    widths = c(4, 2)
  )+  guides(fill = "none")
library(patchwork)

final_fig <- p_full + p_zoom_plot_highlight +
  
  plot_layout(
    widths = c(4, 2),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "top"
  )

final_fig
# PDF (vector; preferred for publication)
ggsave(
  "figures/AF2_tree_zoom_figure.pdf",
  plot = final_fig,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf
)

# TIFF (publication raster)
ggsave(
  "figures/AF2_tree_zoom_figure.tiff",
  plot = final_fig,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# PNG (high-resolution raster)
ggsave(
  "figures/AF2_tree_zoom_figure.png",
  plot = final_fig,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ---------------------------

# Import GADM Nigeria states

# ---------------------------

nigeria_states <- readRDS("data/gadm/gadm41_NGA_1_pk.rds") %>%
  
  st_as_sf()

# ---------------------------

# Select states represented in tree

# ---------------------------

tree_states <- nigeria_states %>%
  
  filter(NAME_1 %in% c(
    
    "Bauchi",
    
    "Kaduna",
    
    "Plateau",
    
    "Abia"
    
  ))

# ---------------------------
# Optional: study site locations
# ---------------------------

nigeria_lgas <- readRDS("data/gadm/gadm41_NGA_2_pk.rds") %>%
  
  st_as_sf()

study_points <- nigeria_lgas %>%
  filter(
    NAME_1 == "Kaduna",
    NAME_2 %in% c("Lere", "Jema'a"))%>%
  st_centroid()

study_coords <- study_points %>%
  mutate(label = NAME_2) %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry()
# ---------------------------
# Map
# ---------------------------

p_tree_map <- ggplot() +
  
  geom_sf(
    data = nigeria_states,
    fill = "grey95",
    colour = "white",
    linewidth = 0.2
  ) +
  geom_sf(
    data = tree_states,
    fill = "#62F962",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_sf_text(
    data = tree_states,
    aes(label = toupper(NAME_1)),
    size = 2.5,
    fontface = "bold"
  ) +
geom_star(
  data = study_coords,
  aes(
    x = X,
    y = Y
  ),
  starshape = 1,
  size = 3,
  fill = "gold",
  colour = "black"
) +
  geom_label(
    data = study_coords,
    aes(
      x = X,
      y = Y,
      label = label
    ),
    nudge_y = 0.25,
    size = 3,
    fontface = "bold.italic",
    fill = "white",
    alpha=0.5,
    colour = "black",
    label.size = 0
  )+
  # geom_sf(
  #   data = study_points,
  #   shape = 8,
  #   size = 3,
  #   fill = "gold",
  #   colour = "black"
  # ) +
  annotation_scale(
        location = "bl",
    width_hint = 0.25
  ) +
  annotation_north_arrow(
    location = "tr",
    style = north_arrow_fancy_orienteering(),
    height = unit(0.5, "cm"),
    width = unit(0.5, "cm")
  ) +
  coord_sf(
    xlim = st_bbox(tree_states)[c("xmin", "xmax")],
    ylim = st_bbox(tree_states)[c("ymin", "ymax")]
  ) +
  theme_void()
# +
#   theme(
#     text = element_text(size = 14)
#   )

p_tree_map


library(patchwork)


p_zoom_plot <- p_zoom_plot +
  theme(
    plot.margin = margin(5, -10, 5, 5)
  )

p_tree_map <- p_tree_map +
  theme(
    plot.margin = margin(5, 5, 5, -10)
  )

final_fig <- p_zoom_plot2+ p_tree_map+ 
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 14,
        face = "bold"
      )
    )
  );final_fig



ggsave(
  "figures/AF2_tree_map.pdf",
  plot = final_fig,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf
)
ggsave(
  "figures/AF2_tree_map.png",
  plot = final_fig,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  "figures/AF2_tree_map.tiff",
  plot = final_fig,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

