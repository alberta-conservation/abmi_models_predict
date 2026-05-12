
library(sf)
library(terra)
library(tidyverse)
library(ggplot2)
library(tidyterra)

# Get the boundary file for the polygon grid used for the predictions
osr <- st_read("0_data/raw/shapefiles/osr_grid_boundary.shp")

# Get the abmi bird list
abmi_birds <- readRDS("0_data/processed/abmiBirdList.rds")

# Get the layer of the individual scheme approvals (lease areas) and filter those without an APP_HOLDER listed
leases <- st_read("0_data/raw/shapefiles/In Situ Oil Sands Scheme Approval_NAD83_10TM_AEPForest.shp") |> 
  mutate(co_prodField = paste(PROD_FIELD, APP_HOLDER, sep = "_"), 
         co_prodPool = paste(PROD_POOL, APP_HOLDER, sep = "_"), 
         area_ha = as.numeric(st_area(geometry))/100^2) |> 
  filter(!is.na(APP_HOLDER) & !is.na(PROD_FIELD))

hist(leases$area_ha[leases$area_ha < 200])

ggplot() + geom_sf(data = osr) + geom_sf(data = leases, aes(fill = PROD_POOL), show.legend = FALSE)


nlevels(as.factor(leases$SCHEME_NO))
leases$area_ha <- as.numeric(st_area(leases))/100^2

nlevels(as.factor(leases$PROD_POOL))

lease100 <- leases |> filter(area_ha >= 100)

lease400 <- leases |> filter(area_ha >= 400)

lease10km <- leases |> filter(area_ha >= 1000)

nlevels(as.factor(lease400$co_prodPool))

tf1 <- table(lease400$co_prodPool[which(lease400$PROD_FIELD == pf[1])])
tf2 <- table(lease400$co_prodPool[which(lease400$PROD_FIELD == pf[2])])
tf3 <- table(lease400$co_prodPool[which(lease400$PROD_FIELD == pf[3])])
tf4 <- table(lease400$co_prodPool[which(lease400$PROD_FIELD == pf[4])])

p <- levels(as.factor(lease400$PROD_POOL))
pf <- levels(as.factor(lease400$PROD_FIELD))
cpf <- levels(as.factor(lease400$co_prodField))
cpp <- levels(as.factor(lease400$co_prodPool))
n <- levels(as.factor(lease400$SCHEME_NO))

ggplot() + geom_sf(data = osr) + geom_sf(data = lease400, aes(fill = co_prodPool), show.legend = FALSE)


## -------------------------------------------------------------------------------------
## Start by including all leases
## -------------------------------------------------------------------------------------

# Create the exposure table for btnw 
btnw_ref <- rast("2_pipeline/store/spp_pred_reference/BlackthroatedGreenWarbler_osr_reference.tif")$Species
plot(btnw_ref)

refpop <- round(sum(values(btnw_ref), na.rm = TRUE))
n <- levels(as.factor(leases$SCHEME_NO))

btnw_exp <- do.call(rbind, lapply(1:length(n), function(x){
  d <- leases |> filter(SCHEME_NO == n[x])
  d1 <- crop(btnw_ref, vect(d), mask = TRUE) # raster clipped to co_prodpool area
  d2 <- round(sum(values(d1), na.rm = TRUE)) # co_prodpool population 
  d3 <- as.numeric(st_area(d))/as.numeric(st_area(osr))*100 # expected percent of osr pop 
  d4 <- d2/refpop*100 # observed percent of osr pop 
  d5 <- (d4/d3 )
  out <- st_as_sf(data.frame(spp = "btnw", osa = d$PROD_FIELD, lease_holder = d$APP_HOLDER, lease = n[x], lease_pop = d2, pop_pct = round(d4, 3), index = round(d5, 3)), geometry = d$geometry)
}))

saveRDS(btnw_exp, "2_pipeline/store/btnw_exposure.rds")
writeVector(vect(btnw_exp), filename = "2_pipeline/store/btnw_exposure.shp")

btnw_current <- rast("2_pipeline/store/spp_pred_current/BlackthroatedGreenWarbler_osr_current.tif")$Species
plot(btnw_current)

btnw_exp_current <- do.call(rbind, lapply(1:length(n), function(x){
  d <- leases |> filter(SCHEME_NO == n[x])
  d1 <- crop(btnw_current, vect(d), mask = TRUE) # raster clipped to co_prodpool area
  d2 <- round(sum(values(d1), na.rm = TRUE)) # co_prodpool population 
  d3 <- as.numeric(st_area(d))/as.numeric(st_area(osr))*100 # expected percent of osr pop 
  d4 <- d2/refpop*100 # observed percent of osr pop 
  d5 <- (d4/d3 )
  out <- st_as_sf(data.frame(spp = "btnw", osa = d$PROD_FIELD, lease_holder = d$APP_HOLDER, lease = n[x], lease_pop = d2, pop_pct = round(d4, 3), index = round(d5, 3)), geometry = d$geometry)
}))

saveRDS(btnw_exp_current, "2_pipeline/store/btnw_exposure_current.rds")


## ------------------------------------------------------------------------------------
## Do the exposure for all species 
## ------------------------------------------------------------------------------------

ref_ls <- list.files("2_pipeline/store/spp_pred_reference", pattern = ".tif")
current_ls <- list.files("2_pipeline/store/spp_pred_current", pattern = ".tif")

osr_exposure_metrics <- function(exp_areas, ref_area, reference = TRUE){
  
  if(reference){
    fol <- "spp_pred_reference"
    ls <- list.files("2_pipeline/store/spp_pred_reference", pattern = ".tif")
    ref = "reference"
  }else{
    fol <- "spp_pred_current"
    ls <- list.files("2_pipeline/store/spp_pred_current", pattern = ".tif")
    ref = "current"
  }
  
  osr_exp <- do.call(rbind, lapply(1:length(ls), function(x){
    sp <- strsplit(ls[x], "_")[[1]][1]
    r <- rast(file.path("2_pipeline/store", fol, ls[x]))$Species
    
    # Calculate the exposure metrics for each lease area 
    refpop <- round(sum(values(r), na.rm = TRUE))
    n <- levels(as.factor(exp_areas$SCHEME_NO))
    
    # Calculate exposure metrics
    spp_exp <- do.call(rbind, lapply(1:length(n), function(x){
      d <- leases |> filter(SCHEME_NO == n[x])
      d1 <- crop(r, vect(d), mask = TRUE) # raster clipped to lease area
      d2 <- round(sum(values(d1$Species), na.rm = TRUE)) # co_prodpool population 
      d3 <- as.numeric(st_area(d))/as.numeric(st_area(ref_area))*100 # expected percent of osr pop 
      d4 <- d2/refpop*100 # observed percent of osr pop 
      d5 <- (d4/d3)
      out <- st_as_sf(data.frame(spp = sp, osa = d$PROD_FIELD, prod_pool = d$PROD_POOL, lease_holder = d$APP_HOLDER, lease = n[x], lease_area_ha = as.numeric(st_area(d))/100^2, osr_pop = refpop, lease_pop = d2, lease_pct = round(d4, 3), lease_index = round(d5, 3)), geometry = d$geometry)
    }))
    
    return(spp_exp)
  }))
  
  return(osr_exp)
}

system.time({
  osr_exp_ref <- osr_exposure_metrics(exp_areas = leases, ref_area = osr, reference = TRUE)
  beep(9)
})

save(osr_exp_ref, file = "2_pipeline/store/abmi_reference_exposure_metrics.rda")

system.time({
  osr_exp_current <- osr_exposure_metrics(exp_areas = leases, ref_area = osr, reference = FALSE)
  beep(9)
})

save(osr_exp_current, file = "2_pipeline/store/abmi_current_exposure_metrics.rda")




