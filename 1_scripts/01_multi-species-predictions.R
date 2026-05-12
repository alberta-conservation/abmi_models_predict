
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


# Remove any grid cells not included in climate data (i.e. too far outside Alberta)
rm <- which(is.na(climate.input$wN))

climate.input.rm <- climate.input[-rm, ]
veg.data.current.rm <- veg.data.current[-which(rownames(veg.data.current) %in% rownames(climate.input)[rm]), ]
veg.data.reference.rm <- veg.data.reference[-which(rownames(veg.data.reference) %in% rownames(climate.input)[rm]), ]





species.lookup <- abmi_species()

birds <- species.lookup |> filter(Taxon == "Birds")
rownames(birds)

saveRDS(birds, "0_data/processed/abmiBirdList.rds")

# spp_list <- c("AlderFlycatcher", "BlackthroatedGreenWarbler", "CanadaWarbler", "Ovenbird", "RubycrownedKinglet", "SwainsonsThrush")
# i <- 0

spp_predict <- function(spp_list, folder, reference = TRUE, i = 0){
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
    model.output <- species_predict(species = spp, 
                                    veg = veg.data, 
                                    climate = climate.input.rm, 
                                    modified = FALSE, 
                                    boot = i)
    if(!is.null(model.output$Vegetation)){
      if(reference){
        truncated.pred <- truncate_predict(spp, current = NULL, reference = model.output$Vegetation)
        abundance <- truncated.pred$reference*10^2
        saveRDS(truncated.pred$reference, paste0("2_pipeline/store/", fol, "/", spp, "_osr_predictions_1km_", ref, ".rds"))
      }else{
        truncated.pred <- truncate_predict(spp, current = model.output$Vegetation, reference = NULL)
        abundance <- truncated.pred$current*10^2
        saveRDS(truncated.pred$current, paste0("2_pipeline/store/", fol, "/", spp, "_osr_predictions_1km_", ref, ".rds"))
      }
      
      # Create raster map from prediction 
      species.pred <- spatial.locations[match(names(abundance), spatial.locations$fid), ]
      species.pred$Species <- abundance
      species.pred <- rast(x = species.pred,
                           type = "xyz")
      
      writeRaster(species.pred, filename = paste0("2_pipeline/store/", fol, "/", spp, "_osr_", ref, ".tif"), overwrite = TRUE)
    }else{
      return(NULL)
    }
    
  })
}

spp_predict(spp_list = birds$SpeciesID[1:2], folder = NULL, reference = TRUE)
spp_predict(spp_list = birds$SpeciesID, folder = NULL, reference = FALSE)

# Get the final list of spp with predictions for the OSR
spp_ls <- list.files("2_pipeline/store/spp_pred_reference", pattern = ".tif")
spp_list <- sapply(1:length(spp_ls), function(x) strsplit(spp_ls[x], "_")[[1]][1])
final_spp_list <- birds |> filter(SpeciesID %in% spp_list)

write.csv(final_spp_list, file = "2_pipeline/store/final_abmi_spp_list.csv")
