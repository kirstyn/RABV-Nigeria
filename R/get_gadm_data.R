# Install packages if needed
# install.packages(c("sf", "ggplot2", "dplyr", "patchwork", "geodata"))

library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(geodata)

# ---------------------------
# Download administrative data
# ---------------------------

# Nigeria states
nigeria_states <- gadm(
  country = "NGA",
  level = 1,
  path = "data"
)
# Store data for future use
st_write(
  nigeria_states,
  "data/spatial/nigeria_states.gpkg",
  delete_dsn = TRUE
)
# Nigeria LGAs
nigeria_lgas <- gadm(
  country = "NGA",
  level = 2,
  path = "data"
)
st_write(
  nigeria_lgas,
  "data/spatial/nigeria_lgas.gpkg",
  delete_dsn = TRUE
)

