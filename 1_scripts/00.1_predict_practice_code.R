
library(veghfsoil)
library(sf)
library(tidyverse)
library(ggplot2)

d <- read.csv("0_data/processed/osr_grid_intersection_bfveg_2022.csv")
g <- st_read("0_data/processed/shapefiles/osr_grid.shp"); gc()


# all.equal(d$RECALCAREA, d$Shape_Area)

## If above result != 'FALSE'
# AREA_COL<-"RECALCAREA"
# if (AREA_COL != "Shape_Area") {
#   d[["Shape_Area"]] <- d[[AREA_COL]]
#   d[[AREA_COL]] <- NULL
# }
# cat("OK\n\n")

# Correct the feature types based on the 2022 HFI classes
d$FEATURE_TY <- as.character(d$FEATURE_TY)
d$FEATURE_TY[d$FEATURE_TY %in% c("TIMBER-HARVEST-GREEN-AREA", "TIMBER-HARVEST-WHITE-AREA")] <- "HARVEST-AREA"
d$FEATURE_TY[d$FEATURE_TY %in% c("OTHER-LINEAR-FEATURES")] <- "TRAIL"
d$FEATURE_TY[d$FEATURE_TY %in% c("REC-TRAIL")] <- "TRAIL"
d$FEATURE_TY[d$FEATURE_TY %in% c("WOODY-VEGETATON-REMOVAL")] <- "WOODY-VEGETATION-REMOVAL"
d$FEATURE_TY <- as.factor(d$FEATURE_TY)

# ggplot() + geom_sf(data = g1) + geom_sf(data = lc1)

# Obtain long summary
d.long <- make_landcover_long(
  landcover = d,
  col.label = "fid_1",
  col.baseyear = 2022,
  col.hfyear = "YEAR",
  col.veg = "Combined_ChgByCWCS",
  col.soil = "Soil_Type_1",
  hf.fine = TRUE,
  burn.cc = TRUE,
  age.correction = "Age_Correction",
  ver.id = "V7.0"
); gc()


# Obtain wide summary
d.wide <- make_landcover_wide(
  long.output = d.long,
  col.label = "fid_1",
  col.area = "Shape_Area",
  hf.fine = TRUE,
  tol = 0,
  sparse = TRUE,
  assign.unknown.ages = TRUE,
  age.data = "Maltman.Old",
  ver.id = "V7.0",
  rm0 = TRUE
); gc()

lc.coef.update <- landcover.coef.lookup

lc.coef.update$Vegetation$COEF <- replace_values(lc.coef.update$Vegetation$COEF, 
                                                 "TreedFenR" ~ "TreedFen", 
                                                 "TreedFen1" ~ "TreedFen", 
                                                 "TreedFen2" ~ "TreedFen", 
                                                 "TreedFen3" ~ "TreedFen", 
                                                 "TreedFen4" ~ "TreedFen", 
                                                 "TreedFen5" ~ "TreedFen", 
                                                 "TreedFen6" ~ "TreedFen", 
                                                 "TreedFen7" ~ "TreedFen", 
                                                 "TreedFen8" ~ "TreedFen", 
                                                 "TreedFen9" ~ "TreedFen", 
                                                 "Urban" ~ "UrbanIndustrial", 
                                                 "Industrial" ~ "IndustrialRural", 
                                                 "Rural" ~ "RuralResidential")


landcover.out <- clean_landcover(
  data.in = as.matrix(d.wide$veg.current),
  landscape.lookup = lc.coef.update,
  type = "Vegetation",
  class.in = "ID",
  class.out = "COEF")

# Convert the areas to proportions
lc.out <- landcover.out/rowSums(landcover.out)
summary(rowSums(lc.out))


library(ABMIexploreR)

# Load the bioclimatic raster into memory
abmi_load_bioclimatic()

cent <- st_centroid(g)

# Extraction of single species coeffcients for vegetation models
single.species <- coefficient_extraction(species = "BlackthroatedGreenWarbler", model = "Vegetation")
str(single.species)

# plot_coef(species = "BlackthroatedGreenWarbler", model = "Vegetation")

spatial.locations <- vect(cent)
crs(spatial.locations) <- "EPSG:3400"

climate.input <- extract_climate(spatial.grid = spatial.locations, 
                                 cell.id = cent$fid,
                                 reproject = FALSE); gc()


## define species and bootstrap id
spp <- "BlackthroatedGreenWarbler"
i <- 0

## Define the compositional data
veg.data <- lc.out

rm <- which(is.na(climate.input$wN))

climate.input.rm <- climate.input[-rm, ]
veg.data.rm <- veg.data[-which(rownames(veg.data) %in% rownames(climate.input)[rm]), ]

save(climate.input.rm, veg.data.rm, spatial.locations, file = "2_pipeline/store/abmi_prediction_data.Rdata")

## Make single species prediction
model.output <- species_predict(species = spp, 
                                veg = veg.data.rm, 
                                climate = climate.input.rm, 
                                modified = FALSE, 
                                boot = i)

saveRDS(model.output, "2_pipeline/store/btnw_osr_predictions_1km_prop.rds")

truncated.pred <- truncate_predict(spp, current = model.output$Vegetation, reference = NULL)
str(truncated.pred)
summary(truncated.pred$current)

## Convert densities to abundance
btnw_abundance <- truncated.pred$current*10^2
summary(btnw_abundance)


# Create raster map from prediction 
species.pred <- spatial.locations[match(names(btnw_abundance), spatial.locations$fid), ]
species.pred$Species <- btnw_abundance
species.pred <- rast(x = species.pred,
                     type = "xyz")
# crs(species.pred) <- "EPSG:3400"

plot_species(spat.raster = species.pred,
             variable = "Species")

writeRaster(species.pred, filename = "2_pipeline/store/spp_pred/btnw_osr.tif")
