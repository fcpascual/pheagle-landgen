pacman::p_load(algatr, geodata, terra, RStoolbox, raster, 
               vcfR, viridis, dplyr, sf, here)

# Download Philippines WorldClim bioclimatic variables
PH_bioclim <- worldclim_country(
  country = "PHL",
  var = "bio",
  path = "worldclim_data"
)

# Download Philippine administrative boundaries
PH_bounds <- gadm(
  country = "PHL",
  level = 1,
  path = "gadm_data"
)

# Group together Mindanao provinces
mindanao_provinces <- c(
  "Agusan del Norte",
  "Agusan del Sur",
  "Bukidnon",
  "Camiguin",
  "Compostela Valley",
  "Davao del Norte",
  "Davao del Sur",
  "Davao Oriental",
  "Lanao del Norte",
  "Lanao del Sur",
  "Maguindanao",
  "Misamis Occidental",
  "Misamis Oriental",
  "North Cotabato",
  "Sarangani",
  "South Cotabato",
  "Sultan Kudarat",
  "Surigao del Norte",
  "Surigao del Sur",
  "Zamboanga del Norte",
  "Zamboanga del Sur",
  "Zamboanga Sibugay"
)

# Select only those provinces
Mindanao <- PH_bounds[
  PH_bounds$NAME_1 %in% mindanao_provinces,
]
Mindanao <- terra::aggregate(Mindanao)
plot(Mindanao)

# Crop bioclim to Mindanao
Mindanao_bioclim <- terra::crop(
  PH_bioclim,
  Mindanao
)

# Mask bioclim to Mindanao shape
Mindanao_bioclim <- terra::mask(
  Mindanao_bioclim,
  Mindanao
)

# Inspect bioclim layer/s on Mindanao map
layer <- 1
plot(Mindanao_bioclim[[layer]])

# Check how many bioclim layers are stored
names(Mindanao_bioclim)
nlyr(Mindanao_bioclim)


# Convert bioclim SpatRaster to RasterStack
Mindanao_bioclim_raster <- raster::stack(
  Mindanao_bioclim
)

# Run rasterPCA
Mindanao_PCA <- RStoolbox::rasterPCA(
  Mindanao_bioclim_raster,
  spca = TRUE
)
Mindanao_PCA # Inspect PCA
summary(Mindanao_PCA$model) # PCA variances
Mindanao_PCA$model$rotation # PCA loadings


# Get first n PCs
nPCs <- 6
Mindanao_env_PCsubset <- Mindanao_PCA$map[[1:nPCs]]

names(Mindanao_env_PCsubset) <- paste0(
  "Mindanao_rPCA",
  seq_len(nPCs)
)
Mindanao_env_PCsubset
plot(Mindanao_env_PCsubset)

# Convert SpatRaster to Rasterstack
# Mindanao_env_rasterstack <- raster::stack(Mindanao_env_PCsubset)
# 
# class(Mindanao_env_rasterstack)
# raster::crs(Mindanao_env_rasterstack)
# raster::extent(Mindanao_env_rasterstack)


# K-means clustering of >3 PCs
pixel_data <- raster::values(Mindanao_PCA$map[[1:nPCs]])
clean_indices <- which(complete.cases(pixel_data))
clean_pixels <- pixel_data[clean_indices, ]

set.seed(327)
km_fit <- kmeans(
  clean_pixels, 
  centers = 10, 
  iter.max = 500, 
  nstart = 5)

Mindanao_env_clusters <- rast(Mindanao_PCA$map[[1]])
values(Mindanao_env_clusters) <- NA

values(Mindanao_env_clusters)[clean_indices] <- km_fit$cluster

plot(Mindanao_env_clusters, 
     type = "classes", 
     col = terra::map.pal("viridis", 10),
     main = "Mindanao Environmental Clusters",
     plg = list(title = "Cluster")
     )


# -----------------------------------------------------
# Load in Philippine Eagle VCF and correct naming
PHEagle_VCF <- read.vcfR("06_BQSR_snps.dataset.recode.vcf")
# colnames(PHEagle_VCF@gt)[2:ncol(PHEagle_VCF@gt)] <- paste0("sample", colnames(PHEagle_VCF@gt)[2:ncol(PHEagle_VCF@gt)])

# Load in sampling coordinates file
PHEagle_coords <- read.csv("coords_long-lat2.csv")

# Convert to sf
PHEagle_coords_sf <- sf::st_as_sf(
  PHEagle_coords,
  coords = c("x", "y"),
  crs = 4326,
  remove = FALSE
)

# Convert coords to planar
## Project sampling coordinates
PHEagle_coords_UTM <- sf::st_transform(
  PHEagle_coords_sf,
  crs = 32651
)

## Project TESS grid
TESS_grid_UTM <- terra::project(
  Mindanao_env_PCsubset[[1]],
  "EPSG:32651"
)


# Check distance between sampling locations
sample_cells <- terra::cellFromXY(
  Mindanao_PCA$map[[1]],
  PHEagle_coords[, c("x", "y")]
)

cell_positions <- data.frame(
  sample = seq_len(nrow(PHEagle_coords)),
  x = PHEagle_coords$x,
  y = PHEagle_coords$y,
  row = terra::rowFromCell(
    Mindanao_PCA$map[[1]],
    sample_cells
  ),
  col = terra::colFromCell(
    Mindanao_PCA$map[[1]],
    sample_cells
  )
)

cell_positions


# --------------------VCF filtering--------------------
# View environmental layers + sampling coordinates
plot(
  Mindanao_env_clusters, 
  col = viridis(100), 
  axes = FALSE
  )

plot(Mindanao_env_clusters, 
     type = "classes", 
     col = terra::map.pal("viridis", 10),
     main = "Mindanao Environmental Clusters & Sampling Locations",
     plg = list(title = "Cluster")
)
sbar(d = 100, xy = "bottomleft" , label = "100 km")
points(PHEagle_coords, pch = 19, col = "#ff1f5b")

# Remove individuals w/o sampling coordinates
cols_to_remove <- c(
  "01", "02", "03", "04", "09", "15", "16", "17", "18", "20",
  "22", "23", "25", "26", "27", "29", "31", "33", "35")
cols_to_keep <- c(
  "FORMAT", 
  "05", "06", "07", "08", "10", "11", "12", "13",
  "14", "19", "21", "24", "28", "30", "32", "34")
PHEagle_VCF_filt <- PHEagle_VCF[, cols_to_keep]

# Remove invariant SNPs after subsetting
poly_vcfr <- vcfR::is.polymorphic(
  PHEagle_VCF_filt,
  na.omit = TRUE
)

# Check how many are polymorphic/invariant
table(poly_vcfr)

# Keep only polymorphic SNPs
PHEagle_VCF_filt <- PHEagle_VCF_filt[poly_vcfr, ]

# Verify filtering
table(
  vcfR::is.polymorphic(
    PHEagle_VCF_filt,
    na.omit = TRUE
  )
)
PHEagle_VCF_filt


# --------------------SAVE TO FILES--------------------
# Filtered master VCF
vcfR::write.vcf(
  PHEagle_VCF_filt,
  file = "ALGATR_input/PHEagle_filt.vcf"
)

# LD pruned VCF
algatr::ld_prune(
  vcf = "ALGATR_input/PHEagle_filt.vcf",
  out_name = "PHEagle_filt",
  out_format = "vcf"
) # --> already writes out a file

# Sampling coordinates
## longitude-latitude
write.csv(
  PHEagle_coords,
  "ALGATR_input/PHEagle_coords_longlat.csv",
  row.names = FALSE
)

## planar/UTM
saveRDS(
  PHEagle_coords_UTM,
  "ALGATR_input/PHEagle_coords_UTM.rds"
)

# Environmental PCs
saveRDS(
  Mindanao_env_PCsubset,
  "ALGATR_input/Mindanao_env_PCsubset.rds"
)

terra::writeRaster(
  Mindanao_env_PCsubset,
  "ALGATR_input/Mindanao_env_PCsubset.tif",
  overwrite = TRUE
)


# Rasters
## Mindanao boundary mask
saveRDS(
  Mindanao,
  "ALGATR_input/Mindanao_boundary.rds"
)

## Mindanao env clusters
saveRDS(
  Mindanao_env_clusters,
  "ALGATR_input/Mindanao_env_clusters.rds"
)
