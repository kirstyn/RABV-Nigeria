# Install packages if needed
# install.packages(c("sf", "ggplot2", "dplyr", "patchwork", "geodata"))

library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(geodata)
library(ggspatial)

# ---------------------------
# Import pre-downloaded gadm data
# ---------------------------
nigeria_states <- readRDS("data/gadm/gadm41_NGA_1_pk.rds")
nigeria_lgas <- readRDS("data/gadm/gadm41_NGA_2_pk.rds")
nigeria_states <- sf::st_as_sf(nigeria_states)
nigeria_lgas  <- sf::st_as_sf(nigeria_lgas)
# ---------------------------
# Select Kaduna State
# ---------------------------

kaduna <- nigeria_states %>%
  filter(NAME_1 == "Kaduna")


# Select Kaduna LGAs
kaduna_lgas <- nigeria_lgas %>%
  filter(NAME_1 == "Kaduna")


# Check LGA names
unique(kaduna_lgas$NAME_2)


# ---------------------------
# Highlight study LGAs
# ---------------------------

study_lgas <- kaduna_lgas %>%
  filter(NAME_2 %in% c("Lere", "Jema'a"))


# ---------------------------
# Panel A: Nigeria map
# ---------------------------
kaduna_bbox <- st_as_sfc(st_bbox(kaduna))
p1 <- ggplot() +
  geom_sf(data = nigeria_states,
          fill = "grey90",
          colour = "white") +
  geom_sf(data = kaduna,
          fill = "#D6604D",
          colour = "black") +
  geom_sf(data = kaduna_bbox,
          fill = NA,
          colour = "black",
          linewidth = 0.8) +
  geom_sf_text(data = nigeria_states %>% filter(NAME_1 != "Kaduna"),
               aes(label = NAME_1),
               size = 4,
               check_overlap = TRUE) +
  geom_sf_text(data = kaduna,
               aes(label = NAME_1),
               size = 6,
               fontface = "bold") +
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
  theme_minimal() +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, colour = "grey40"),
    axis.title = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    panel.background = element_blank()
  ) +
  ggtitle("A");p1


# ---------------------------
# Panel B: Kaduna LGAs
# ---------------------------

p2 <- ggplot() +
  geom_sf(data = kaduna_lgas,
          fill = "grey95",
          colour = "white") +
  geom_sf(data = study_lgas,
          fill = "#D6604D",
          colour = "black") +
  geom_sf_text(data = kaduna_lgas %>% 
                 filter(!NAME_2 %in% c("Lere", "Jema'a")),
               aes(label = NAME_2),
               size = 5,
               check_overlap = TRUE) +
  geom_sf_text(data = study_lgas,
               aes(label = NAME_2),
               size = 6,
               fontface = "bold") +
  annotation_scale(
    location = "bl",
    width_hint = 0.25
  ) +
  annotation_north_arrow(
    location = "tr",
    style = north_arrow_fancy_orienteering(),
    height = unit(1, "cm"),
    width = unit(1, "cm")
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
  )+
  ggtitle("B");p2



# ---------------------------
# Combine
# ---------------------------

figure1 <- p1 +p2 
figure1


ggsave(
  "figures/Fig1_Kaduna_rabies_study_map.png",
  figure1,
  width = 12,
  height = 5,
  dpi = 600
)
ggsave(
  "figures/Fig1_Kaduna_rabies_study_map.tiff",
  figure1,
  width = 12,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

