
library(veghfsoil)
library(sf)
library(tidyverse)
library(ggplot2)

d <- st_read("0_data/processed/g2_int_test.gdb", "g2_int_test")
g <- st_read("0_data/processed/g2_grid_test.shp")

save(d, g, file = "0_data/processed/g2_test_data.Rdata")

a <- as.numeric(st_area(d))

all.equal(a, d$Shape_Area)

d$RECALCAREA <- as.numeric(st_area(d))

AREA_COL<-"RECALCAREA"
if (AREA_COL != "Shape_Area") {
  d[["Shape_Area"]] <- d[[AREA_COL]]
  d[[AREA_COL]] <- NULL
}
cat("OK\n\n")

lc <- d

# Correct the feature types based on the 2022 HFI classes
lc$FEATURE_TY <- as.character(lc$FEATURE_TY)
lc$FEATURE_TY[lc$FEATURE_TY %in% c("TIMBER-HARVEST-GREEN-AREA", "TIMBER-HARVEST-WHITE-AREA")] <- "HARVEST-AREA"
lc$FEATURE_TY[lc$FEATURE_TY %in% c("OTHER-LINEAR-FEATURES")] <- "TRAIL"
lc$FEATURE_TY[lc$FEATURE_TY %in% c("REC-TRAIL")] <- "TRAIL"
lc$FEATURE_TY[lc$FEATURE_TY %in% c("WOODY-VEGETATON-REMOVAL")] <- "WOODY-VEGETATION-REMOVAL"
lc$FEATURE_TY <- as.factor(lc$FEATURE_TY)

lc1 <- lc |> filter(fid_1 %in% as.numeric(levels(as.factor(d$fid_1))[1:100]))

g1 <- g |> filter(fid %in% lc1$fid_1)

ggplot() + geom_sf(data = g1) + geom_sf(data = lc1)

lc2 <- st_drop_geometry(lc1)

# Obtain long summary
d.long <- make_landcover_long(
  landcover = lc2,
  col.label = "fid_1",
  col.baseyear = 2022,
  col.hfyear = "YEAR",
  col.veg = "Combined_ChgByCWCS",
  col.soil = "Soil_Type_1",
  hf.fine = TRUE,
  burn.cc = TRUE,
  age.correction = "Age_Correction",
  ver.id = "V7.0"
)


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
)

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

library(ABMIexploreR)

# Load the bioclimatic raster into memory
abmi_load_bioclimatic()

cent <- st_centroid(g1)

# Extraction of single species coeffcients for vegetation models
single.species <- coefficient_extraction(species = "BlackthroatedGreenWarbler", model = "Vegetation")
str(single.species)

plot_coef(species = "BlackthroatedGreenWarbler", model = "Vegetation")

spatial.locations <- vect(st_coordinates(cent), type = "point")
crs(spatial.locations) <- "EPSG:3400"

climate.input <- extract_climate(spatial.grid = spatial.locations, 
                                 cell.id = cent$fid,
                                 reproject = FALSE)


## define species and bootstrap id
spp <- "BlackthroatedGreenWarbler"
i <- 0

## Define the compositional data
veg.data <- landcover.out

rm <- which(is.na(climate.input$wN))

climate.input.rm <- climate.input[-rm, ]
veg.data.rm <- veg.data[-which(rownames(veg.data) %in% rownames(climate.input)[rm]), ]

## Make single species prediction
model.output <- species_predict(species = spp, 
                                veg = veg.data.rm, 
                                climate = climate.input.rm, 
                                modified = FALSE, 
                                boot = i)

saveRDS(model.output, "2_pipeline/store/btnw_g2_test_predictions.rds")



