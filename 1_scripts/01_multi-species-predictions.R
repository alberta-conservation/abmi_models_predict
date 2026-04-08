
library(sf)
library(terra)
library(tidyverse)
library(ABMIexploreR)
library(veghfsoil)


load("2_pipeline/store/abmi_prediction_data.Rdata")

# Format the data
spatial.locations <- vect(cent)

## Define the compositional data
veg.data.current <- lc.out.current
veg.data.reference <- lc.out.reference


# Remove any grid cells no included in climate data (i.e. too far outside Alberta)
rm <- which(is.na(climate.input$wN))

climate.input.rm <- climate.input[-rm, ]
veg.data.current.rm <- veg.data.current[-which(rownames(veg.data.current) %in% rownames(climate.input)[rm]), ]
veg.data.reference.rm <- veg.data.reference[-which(rownames(veg.data.reference) %in% rownames(climate.input)[rm]), ]





species.lookup <- abmi_species()

birds <- species.lookup |> filter(Taxon == "Birds")
rownames(birds)

spp_list <- c("AlderFlycatcher", "BlackthroatedGreenWarbler", "CanadaWarbler", "Ovenbird", "RubycrownedKinglet", "SwainsonsThrush")
i <- 0


spp_predict <- function(spp_list, folder, reference){
  nspec <- length(spp_list)
  
  if(reference){
    fol <- "spp_pred_reference"
    veg.data <- veg.data.reference.rm
    ref = "reference"
  }else{
    fol <- "spp_pred_current"
    veg.data <- veg.data.current.rm
    ref = "current"
  }
  
  lapply(1:nspec, function(x){
    spp <- spp_list[x]
    # single.species <- coefficient_extraction(species = spp, model = "Vegetation")
    model.output.current <- species_predict(species = spp, 
                                    veg = veg.data, 
                                    climate = climate.input.rm, 
                                    modified = FALSE, 
                                    boot = i)
    if(reference){
      truncated.pred <- truncate_predict(spp, current = NULL, reference = model.output$Vegetation)
      abundance <- truncated.pred$reference*10^2
    }else{
      truncated.pred <- truncate_predict(spp, current = model.output$Vegetation, reference = NULL)
      abundance <- truncated.pred$current*10^2
    }
    
    saveRDS(truncated.pred$current, paste0("2_pipeline/store/", fol, "/", spp, "_osr_predictions_1km_", ref, ".rds"))
    
    
    # Create raster map from prediction 
    species.pred <- spatial.locations[match(names(abundance), spatial.locations$fid), ]
    species.pred$Species <- abundance
    species.pred <- rast(x = species.pred,
                         type = "xyz")
    
    plot_species(spat.raster = species.pred,
                 variable = "Species")
    
    writeRaster(species.pred, filename = paste0("2_pipeline/store/", fol, "/", spp, "_osr_", ref, ".tif"), overwrite = TRUE)
  })
}

spp_predict(spp_list = spp_list, folder = NULL, reference = TRUE)
spp_predict(spp_list = spp_list, folder = NULL, reference = FALSE)

lapply(1:nspec, function(x){
  spp <- spp_list[x]
  # single.species <- coefficient_extraction(species = spp, model = "Vegetation")
  model.output <- species_predict(species = spp, 
                                  veg = veg.data.current.rm, 
                                  climate = climate.input.rm, 
                                  modified = FALSE, 
                                  boot = i)
  
  truncated.pred <- truncate_predict(spp, current = model.output$Vegetation, reference = NULL)
  saveRDS(truncated.pred$current, paste0("2_pipeline/store/spp_pred_current/", spp, "_osr_predictions_1km_current.rds"))
  
  abundance <- truncated.pred$current*10^2
  
  # Create raster map from prediction 
  species.pred <- spatial.locations[match(names(abundance), spatial.locations$fid), ]
  species.pred$Species <- abundance
  species.pred <- rast(x = species.pred,
                       type = "xyz")

  plot_species(spat.raster = species.pred,
               variable = "Species")
  
  writeRaster(species.pred, filename = paste0("2_pipeline/store/spp_pred_current/", spp, "_osr_current.tif"))
})

## define species and bootstrap id
spp <- "BlackthroatedGreenWarbler"
i <- 0

# Extraction of single species coeffcients for vegetation models
single.species <- coefficient_extraction(species = spp, model = "Vegetation")
str(single.species)


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





